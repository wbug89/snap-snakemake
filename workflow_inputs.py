from __future__ import annotations

import os
import re
from pathlib import Path
from collections import defaultdict


def _abs(path):
    return os.path.abspath(os.path.expanduser(os.path.expandvars(str(path))))


def safe_id(value):
    value = str(value)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", value):
        raise ValueError(
            f"Unsafe sample/unit ID '{value}'. IDs used in paths must start with an alphanumeric "
            "character and contain only letters, numbers, '.', '_' and '-'."
        )
    return value


def directories(config, key):
    vals = config.get(key, []) or []
    if isinstance(vals, (str, os.PathLike)):
        vals = [vals]
    out = []
    for p in vals:
        ap = _abs(p)
        if not os.path.isdir(ap):
            raise ValueError(f"Configured {key} directory does not exist or is not a directory: {ap}")
        out.append(ap)
    if not out:
        raise ValueError(f"Config '{key}' must contain at least one directory for this stage.")
    return out


def _dict_path(reference):
    p = Path(reference)
    if p.suffix in {'.fa', '.fasta', '.fna'}:
        return str(p.with_suffix('.dict'))
    return str(p) + '.dict'


def reference_inputs(reference, stage):
    reference = _abs(reference)
    fai = reference + '.fai'
    dct = _dict_path(reference)
    required = [reference, fai, dct]
    if stage in {'preprocess', 'full'}:
        required += [reference + x for x in ('.amb', '.ann', '.bwt', '.pac', '.sa')]
    missing = [x for x in required if not os.path.isfile(x) or os.path.getsize(x) == 0]
    if missing:
        raise ValueError('Missing/empty reference prerequisite(s): ' + ', '.join(missing))

    contigs = []
    with open(fai) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                raise ValueError(f"Malformed FASTA index line in {fai}: {line.rstrip()}")
            contigs.append((parts[0], int(parts[1])))
    if not contigs:
        raise ValueError(f"No contigs found in FASTA index: {fai}")
    return required, contigs


def _read_id_file(path):
    if not path:
        return set()
    ap = _abs(path)
    if not os.path.isfile(ap):
        raise ValueError(f"Sample list file does not exist: {ap}")
    out = set()
    with open(ap) as fh:
        for raw in fh:
            s = raw.strip()
            if not s or s.startswith('#'):
                continue
            out.add(s.split()[0])
    return out


def selected_samples(config, available):
    available = {str(x) for x in available}
    include = {str(x) for x in (config.get('include_biological_samples', []) or [])}
    include |= _read_id_file(config.get('include_biological_samples_file'))
    exclude = {str(x) for x in (config.get('exclude_biological_samples', []) or [])}
    exclude |= _read_id_file(config.get('exclude_biological_samples_file'))

    chosen = set(available if not include else available & include)
    chosen -= exclude
    unknown_includes = include - available
    if unknown_includes:
        raise ValueError(
            'Requested biological samples were not found in the current stage inputs: ' +
            ', '.join(sorted(unknown_includes))
        )
    return sorted(chosen)


def _read_assignment(path):
    ap = _abs(path)
    if not os.path.isfile(ap):
        raise ValueError(f"sample_assignment_file does not exist: {ap}")
    assignment = {}
    with open(ap) as fh:
        header_seen = False
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = re.split(r'\t+|\s+', line)
            if not header_seen and parts[0].lower() in {'biological_sample', 'sample'}:
                header_seen = True
                continue
            header_seen = True
            if len(parts) < 2:
                raise ValueError(f"Malformed assignment line {lineno} in {ap}: {line}")
            sample, partition = parts[0], parts[1]
            if sample in assignment:
                raise ValueError(f"Duplicate biological sample in assignment TSV: {sample}")
            assignment[sample] = partition
    return assignment


def discover_marked(config, get_sm_name, reference_contigs=None):
    dirs = directories(config, 'marked_dirs')
    bams_by_sample = defaultdict(list)
    bais_by_sample = defaultdict(list)
    seen_paths = set()

    for d in dirs:
        for p in sorted(Path(d).glob('*.marked.bam')):
            ap = str(p.resolve())
            if ap in seen_paths:
                continue
            seen_paths.add(ap)
            stem = p.name[:-len('.marked.bam')]
            if '__UNIT__' in stem:
                sm = stem.split('__UNIT__', 1)[0]
            else:
                sm = get_sm_name(stem)
            safe_id(sm)
            bai1 = ap + '.bai'
            bai2 = ap[:-4] + '.bai'
            bai = bai1 if os.path.isfile(bai1) else bai2 if os.path.isfile(bai2) else None
            if bai is None:
                raise ValueError(f"Missing BAM index for {ap}; expected {bai1} or {bai2}")
            bams_by_sample[sm].append(ap)
            bais_by_sample[sm].append(bai)

    if not bams_by_sample:
        raise ValueError('No *.marked.bam files were found in marked_dirs.')

    available = set(bams_by_sample)
    chosen = set(selected_samples(config, available))

    assign_path = config.get('sample_assignment_file')
    partition = str(config.get('calling_partition') or '').strip()
    if bool(assign_path) != bool(partition):
        raise ValueError('sample_assignment_file and calling_partition must be set together.')
    if assign_path:
        assignment = _read_assignment(assign_path)
        unassigned = available - set(assignment)
        if unassigned:
            raise ValueError('Discovered samples missing from sample_assignment_file: ' + ', '.join(sorted(unassigned)))
        requested = {s for s, part in assignment.items() if part == partition}
        missing_inputs = requested - available
        if missing_inputs:
            raise ValueError(
                f"Samples assigned to partition {partition} have no marked BAM input: " +
                ', '.join(sorted(missing_inputs))
            )
        chosen &= requested

    if not chosen:
        raise ValueError('No biological samples remain after gVCF-stage selection/assignment.')

    return (
        {sm: sorted(bams_by_sample[sm]) for sm in sorted(chosen)},
        {sm: sorted(bais_by_sample[sm]) for sm in sorted(chosen)},
    )


def discover_gvcfs(config, reference_contigs=None):
    dirs = directories(config, 'gvcf_dirs')
    by_sample = {}
    for d in dirs:
        candidates = list(Path(d).glob('*.g.vcf.gz')) + list(Path(d).glob('*.gvcf.gz'))
        for p in sorted(candidates):
            name = p.name
            if name.endswith('.g.vcf.gz'):
                sm = name[:-len('.g.vcf.gz')]
            else:
                sm = name[:-len('.gvcf.gz')]
            safe_id(sm)
            ap = str(p.resolve())
            if sm in by_sample:
                raise ValueError(
                    f"Duplicate biological sample gVCF '{sm}' found in multiple gvcf_dirs. "
                    "Technical replicates must be merged before HaplotypeCaller, not during joint genotyping."
                )
            if not os.path.isfile(ap + '.tbi'):
                raise ValueError(f"Missing tabix index for gVCF: {ap}.tbi")
            by_sample[sm] = ap

    if not by_sample:
        raise ValueError('No *.g.vcf.gz / *.gvcf.gz files were found in gvcf_dirs.')
    chosen = selected_samples(config, set(by_sample))
    if not chosen:
        raise ValueError('No biological samples remain after joint-stage selection.')
    return {sm: by_sample[sm] for sm in chosen}
