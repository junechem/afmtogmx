# TASKS.md Archive

This directory contains archived TASKS.md files from completed work periods.

## Archive Naming Convention

Filename format: `TASKS_YYYY-MM-DD_to_YYYY-MM-DD.md`

Example: `TASKS_2025-12-03_to_2025-12-05.md`

## Archive File Format

Each archived file must include a greppable summary line at the very top:

```markdown
<!-- ARCHIVE_SUMMARY: Brief description of work done (YYYY-MM-DD to YYYY-MM-DD) -->
```

## Quick Search

To see all archived task summaries:

```bash
grep "ARCHIVE_SUMMARY:" archive/*.md
```

## How to Archive TASKS.md

1. **Determine date range**:
   - Start date: When current TASKS.md started (first task date)
   - End date: Today's date or last task completion date

2. **Create archive file**:
   ```bash
   cp TASKS.md archive/TASKS_START-DATE_to_END-DATE.md
   ```

3. **Add summary line** to the archived file at the very top:
   ```markdown
   <!-- ARCHIVE_SUMMARY: Config system, regression tests, charge removal (2025-12-03 to 2025-12-05) -->
   ```

4. **Create fresh TASKS.md**:
   ```bash
   cat > TASKS.md << 'EOF'
   # TASKS.md

   ## In Progress

   (No active tasks)

   ## Future Work

   (Add future tasks here)
   EOF
   ```

5. **Commit to git**:
   ```bash
   git add archive/ TASKS.md
   git commit -m "Archive completed tasks from START-DATE to END-DATE"
   git push
   ```

## Example Archive File

See `TASKS_2025-12-03_to_2025-12-05.md` for an example of the proper format.
