# Managing SNaP with Git on macOS

The recommended model is:

- **Mac**: primary development and Git/GitHub management.
- **GitHub (`snap-snakemake`)**: central remote repository.
- **PC1 / PC3**: analysis hosts that pull tagged or approved workflow versions.

## First local repository

```bash
cd ~/Projects/snap-snakemake
git init
git branch -M main
git status
```

Review `.gitignore` before staging files. Local run-specific configs, credentials, logs, `.snakemake/`, FASTQ/BAM/VCF files, and workflow outputs should not be committed.

```bash
git add .
git status
git commit -m "Initial release: SNaP v1.0.0"
git tag -a v1.0.0 -m "SNaP v1.0.0"
```

After creating an empty GitHub repository named `snap-snakemake`, add the remote and push:

```bash
git remote add origin git@github.com:<YOUR_GITHUB_USERNAME>/snap-snakemake.git
git push -u origin main
git push origin v1.0.0
```

## Daily development

```bash
git status
git diff
git add <files>
git commit -m "Describe the change"
git push
```

For a new release, update `VERSION` and `CHANGELOG.md`, commit the release, then create an annotated tag.
