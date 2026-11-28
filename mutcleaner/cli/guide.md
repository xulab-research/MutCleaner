# MutCleaner Guide

## User guide

Use `inspect` to examine the input and `clean` to write the results. After cleaning, use `inspect --verify` to check the standard output. Run `mutcleaner guide TOPIC` to read a specific section. If you are an agent handling a cleaning request, start by creating the internal checklist in the Agent workflow below.

From the MutCleaner source directory, you can also use `python -m mutcleaner.cli` in place of `mutcleaner`.

<!-- topic: inspect -->
### Inspect the data

```text
mutcleaner inspect input_directory
mutcleaner inspect references.fasta
mutcleaner inspect input.csv
mutcleaner inspect input.csv --sequence-type protein --reference-id-column protein_id --reference-fasta references.fasta --mutation-column change --variant-sequence-column variant --label-column ddG --condition-column temperature
```

Basic inspection lists directory contents, FASTA IDs and sequence lengths, or table columns, row counts, column statistics, and examples. Table statistics cover the entire table by default. For a quick check, use `--profile-rows N` to profile only the first N rows. Check `profiled_rows` and `full_table_profiled` to verify the scope. `--sample-rows 0` hides sample rows while preserving the profiling scope. Cleaning processes the entire input.

For Excel workbooks with multiple sheets, inspection lists the sheets if none is selected. Use `--sheet` or `--sheet-index` to inspect each sheet. Specify `--sequence-type` explicitly when checking candidate sequences or mutation annotations.

Specify candidate columns to check ID matches, conflicting sequences under the same reference ID, sequence lengths and characters, numeric measurements, and duplicate groups. With `--mutation-column`, inspection evaluates offsets 0 and 1 by default. Repeat `--offset` to test other values, or use the same `--offset-for`, `--position-map`, and `--reference-slice` options as in cleaning. Compare counts of valid positions, matching reference residues, and matching variant sequences for the entire table and for each reference. The report shows details for up to 50 reference groups and flags any truncation. Inspection uses the same default alphabet as cleaning. Enable special residues explicitly when appropriate. If you have annotations and variant sequences only, you can also check the number of reference inference conflicts under each candidate numbering scheme. Interpret each group's counts in light of the data source. During annotation inspection, omitting a wild-type (WT) option makes `WT` a candidate marker. The report distinguishes a candidate marker, an explicit marker, and `--no-wild-type-marker`.

With `--mutation-column`, `inspect` looks for separators in observed multi-substitution annotations before checking individual rows. `annotations.mutation_separator` reports the selected separator, record counts for each candidate, and examples. If you specify `--mutation-separator`, inspection validates annotations against that value. If it finds multiple valid separators and you have specified none, it reports `mixed` and skips row-level validation. Review the formats before normalizing them or inspecting them separately. If every annotation describes a single substitution, continue with validation. Pass the verified separator explicitly to `clean`.

<!-- topic: clean -->
### Choose a cleaning method

```bash
mutcleaner clean input.csv output \
  --sequence-type protein \
  --reference-id-column protein_id \
  --reference-fasta references.fasta \
  --mutation-column change --mutation-separator ',' --offset 1 \
  --wild-type-marker WT \
  --label-column ddG
```

| Available inputs | How to run the CLI |
| --- | --- |
| Reference and substitution annotations | Provide the reference source, `--mutation-column`, and the numbering rules to reconstruct variant sequences. |
| Reference and variant sequences | Provide the reference source and `--variant-sequence-column` to compare equal-length sequences. |
| Substitution annotations and variant sequences | Provide `--infer-reference`, `--mutation-column`, `--variant-sequence-column`, and the numbering rules to infer a reference for each ID. |
| All three | Provide all three inputs to reconstruct variant sequences and compare them with the observed sequences. |

When using mutation annotations, including those built with `--mutation-parts`, explicitly specify the mutation separator and WT marker policy.

Use `--reference-id-column` to identify the reference ID column or `--reference-id` to name a single-reference dataset. To extract an ID from a composite value, supply a single capture group with `--reference-id-pattern`. Cleaning requires exactly one reference source. Use a table column, a FASTA file, a CSV, TSV, or JSON reference table, an explicit sequence, `--reference-from-wt`, or `--infer-reference`. With mutation annotations, pass `--wild-type-marker MARKER` or `--no-wild-type-marker`. Marker matching ignores case and surrounding whitespace. `--reference-from-wt` requires a matching marker in the input. `--infer-reference` requires annotations and observed variant sequences. By default, FASTA IDs come from the first word of each header. Use `--reference-file-id-pattern` to extract IDs from an external reference file. For `--reference-table`, specify the ID and sequence columns for CSV or TSV, or supply a JSON mapping from IDs to sequences.

To infer substitutions from sequences alone, explicitly provide the reference source and `--variant-sequence-column`. When using substitution annotations, provide the available sequence columns. The offset is the number assigned to the first position of the reference in the original numbering scheme. Use `--offset-for ID=OFFSET` for a specific reference. For nonlinear numbering, use `--position-map` with a table containing `reference_id`, `input_position`, and zero-based `reference_index`. The map takes precedence. Unmapped positions use the reference-specific offset or the global offset. If any valid annotations remaining after the initial filters lack a numbering rule, cleaning reports the reference ID and the uncovered positions. `--reference-slice START STOP` takes a zero-based, half-open slice. After slicing, annotations and variant sequences are interpreted relative to that slice.

**Files and columns.** The CLI supports CSV, TSV, TXT, Excel (`.xlsx` and `.xls`), and Parquet. Use `--delimiter` and `--encoding` for text files as needed. TXT defaults to tab separation. For an Excel workbook with multiple sheets, select one with `--sheet` or `--sheet-index` before cleaning. IDs stored as strings retain leading zeros. Convert numeric cells to text first if their leading zeros come from display formatting.

**Table preparation.** The CLI applies `--rename OLD=NEW`, `--constant COLUMN=VALUE`, `--join-column TARGET SOURCES SEPARATOR`, `--split-column SOURCE SEPARATOR TARGETS`, `--replace COLUMN PATTERN REPLACEMENT`, and `--extract TARGET SOURCE PATTERN` in that order. Each option can be repeated. Separate the column names in `SOURCES` or `TARGETS` with commas. `--mutation-parts WT POSITION MUTANT` combines three columns into a substitution annotation. `--filter 'temperature>=25'` performs a numeric comparison, while `=` and `!=` compare text values. Use `--exclude COLUMN=PATTERN` to exclude values matching a regular expression and `--drop-missing COLUMN` to drop rows missing a value in that column. These filters run before reference and WT matching.

**Mutations and sequences.** Use `--sequence-type protein` for amino acid substitutions or `dna` or `rna` for sequences of complete codons. Positions are in amino acids or codons, respectively. Multiple substitutions in a variant stay in the same measurement record. Set `--mutation-separator` to the verified separator. If you have confirmed that every annotation describes a single substitution, you can pass a comma explicitly. If separators are mixed, examine the formats before normalizing or handling them separately. Protein options control whether to allow ambiguous residues and U/O. Protein validation rejects `*`. Use `--min-substitutions`, `--max-substitutions`, and `--mutation-effect` to filter by substitution count and effect. Codon output can represent stop codons. Use `--filter-stop` to exclude records that introduce a stop codon in DNA or RNA. This workflow handles substitutions and equal-length sequence comparisons. Insertions, deletions, and frameshifts require separate handling. Use `--workers` to set concurrency for parallel sequence processing. The default is 4. Use `--to-protein` to translate codon data for protein output. Protein validation rejects translated sequences containing `*`, and variants with no amino acid changes are reported as filtered. When different codon variants translate to the same protein variant, `--aggregate` determines how to handle their measurements.

**Custom table processing.** If the standard options cannot express a dataset-specific operation, pass the same `--custom-cleaner custom.py` to `inspect` and `clean`. The script runs on the machine executing MutCleaner. Run only trusted, reviewed scripts because this option can execute arbitrary Python code.

Define a top-level `custom_cleaner(df)` function. It receives the full table and returns a pandas DataFrame. It can add or remove columns, change values, and filter, reorder, or expand rows. It runs before standard column processing, reference matching, mutation validation, and label processing. Select a sheet first for workbooks with multiple sheets. This example converts an encoding such as `sub(A,12,V)` into a substitution annotation.

```python
def custom_cleaner(df):
    return df.assign(mutation=df["encoded_variant"].str.replace(r"^sub\(([A-Z]),(\d+),([A-Z])\)$", r"\1\2\3", regex=True))
```

```bash
mutcleaner inspect raw.csv --custom-cleaner custom.py
mutcleaner clean raw.csv output --custom-cleaner custom.py \
  --sequence-type protein --reference-id-column protein_id \
  --reference-fasta references.fasta --mutation-column mutation \
  --mutation-separator ',' --no-wild-type-marker --offset 1 --label-column ddG
```

Keep the function dependent only on the input table and packages available in the environment. Avoid external files, network access, and random state. `--profile-rows` limits only the statistics collected after custom processing. The JSON `custom` section reports the SHA-256 of the exact source executed, row and column counts before and after processing, and columns added or removed. Values in retained columns may also change, and the row count difference shows only the net change. `source_row` is the one-based row number in the table after custom and standard table processing.

<!-- topic: measurements -->
### Handle experimental measurements

- Repeat `--label-column` to select multiple metrics. Each metric is exported separately and uses only finite numeric values.
- `--aggregate mean` averages comparable repeated measurements for the same reference. When duplicate groups exist, specify `--aggregate` explicitly. Use `keep` to preserve individual measurements or `first` to take the first record in source order. To average repeats, use `--group-by` or `--split-by` to distinguish experimental conditions, or use `--assume-comparable` after verifying that the measurements are comparable. `first` also requires grouping or confirmation of comparability. If `--group-by` alone would place measurements of the same variant under different conditions in the same standard table, the CLI will ask you to use `--split-by` or filter first.
- Repeat `--split-by` to write a separate dataset for each combination of specified conditions. Both grouping and splitting use the original text values of the condition columns.
- `--aggregate nearest --nearest-by temperature=25 --nearest-by pH=7` compares temperature distance first, then uses pH to break ties. `--nearest-weight` multiplies each distance by its weight within the existing priority order. Unspecified weights default to 1. Positive weights generally leave the priority order unchanged. Missing conditions are treated as infinitely distant for that criterion, and earlier criteria still take precedence. If all distances tie, the first record in source order is kept. Use this strategy only when the conditions and measurement methods are comparable.
- `--subtract-wt` matches actual WT measurements by reference ID and grouping fields. It subtracts the matched WT baseline from variant measurements before aggregating them. Fields used for nearest selection also participate in WT matching. Add baseline grouping fields with `--wt-group-by batch`. If the matching group has one finite WT value, the CLI uses it directly. If it has several, choose `--wt-baseline mean` to average them or `--wt-baseline unique` to exclude records without a unique baseline. When you have only sequence pairs and have confirmed that identical sequences represent measured WT records, pass `--wt-from-equal`. An inferred reference cannot serve as a measured baseline.

Combining measurements across temperatures, pH values, substrates, or batches can change the meaning of the label. If the conditions are unclear, investigate, group by condition, or ask the user. The report accounts for missing WT baselines and unavailable nearest-condition measurements in their respective categories.

<!-- topic: output -->
### Verify the results

With one metric and no splits, results go to `output/<reference_id>/`. Each directory contains `data.csv`, `wt.fasta`, `metadata.json`, and a generated `cleaner.py`. If the underlying export fails for some records, it uses the existing `failed_data.csv` format. For multiple metrics or condition splits, see the script and the JSON report for the mapping from directory numbers to labels and condition values.

Each copy of `cleaner.py` reproduces the entire cleaning run and generates all its outputs. Pass the original input and a new output directory. Use the same MutCleaner build and make sure `mutcleaner` is on your PATH. Edit the path variables in the script if an external reference file or position map has moved.

If you used `--custom-cleaner`, `cleaner.py` embeds the exact source that ran. Rerunning the generated script restores the source to a temporary file and uses it again.

```text
python output/<reference_id>/cleaner.py raw.csv new_output
```

The generated script records the parameters, the task associated with its directory, run statistics, and a verification summary. See `metadata.json` for coverage and output row counts and the JSON report for full summaries by task. If no valid variants remain, the command reports this in JSON, exits with a nonzero status, and leaves the output path uncreated.

```text
mutcleaner inspect output --verify
```

`--save-format dataframe` and `--save-format pickle` use native MutCleaner formats with the base name `dataset`. Each reference directory still includes `cleaner.py`. The verification command above applies to the default standard output.

With multiple labels or condition splits, you can verify the output root, a task directory, or an individual reference directory. The CLI reports missing standard files and checks substitution annotations, references, variant sequences, coverage, and position units. Also check that the selected labels, units, transformations, and repeat handling match the experimental question. Cleaning logs go to stderr, and the JSON report goes to stdout.

<!-- topic: agent -->
## Agent workflow

Cleaning usually means organizing, validating, and exporting existing measurements. Start by checking the columns, references, and meaning of the measurements against the source table and supporting material. **You must confirm with the user how to handle unusual cases. Do not make those decisions silently or treat an explicit CLI option as user approval.** For unresolved rules that affect the result, use the evidence to present viable choices, your recommendation, and the reason for it. Wait for confirmation before applying those rules. You can set technical parameters such as file encoding, separators, and worker count after checking what is appropriate.

**Mandatory fact verification**

Before presenting any substantive factual claim, finding, recommendation, or decision question to the user, complete the following verification cycle internally:

1. **Draft.** Form a preliminary conclusion without presenting it to the user.
2. **Plan.** Break the conclusion into individual factual claims and formulate neutral, independently answerable verification questions. Do not embed the expected answers in these questions.
3. **Verify.** Answer each question independently using actual source material, raw data, tool results, or executable checks. Do not use the preliminary conclusion, model knowledge, or another unverified answer as evidence. Actively check for contradictory evidence and exceptions. Match the scope of verification to the scope of the claim.
4. **Revise.** Compare the verification results with the preliminary conclusion. Correct, narrow, or remove unsupported claims. If evidence remains insufficient, investigate further or explicitly report the uncertainty.

Do not replace verification with reflection on your own previous answer. A plausible explanation is not evidence. Never claim that a source was read or a check was performed unless it actually was.

Distinguish verified source facts, observed data facts, inferences, and unresolved questions. Provide traceable evidence for factual claims that affect cleaning decisions.

Repeat this cycle whenever new facts are introduced or previous conclusions change. Prefer additional investigation and tool calls over unsupported assumptions. Complete verification before presenting conclusions, rather than after the user challenges them.

This rule applies throughout the workflow, including intermediate reports, confirmation questions, execution plans, and final reports. Keep preliminary conclusions, verification questions, and internal checklist updates private; present the verified findings, relevant evidence, and unresolved limitations to the user.

When independent subagents are available and their use is appropriate, you may delegate verification questions without providing the preliminary conclusion or expected answers. Provide only the neutral question and the relevant source material or data access. Subagents are optional; when unavailable, perform the checks yourself.

Before inspecting the input or reading the sources, **create an internal checklist using the template below**. Use an internal task tool when available. Otherwise, track the same items privately as you work. Keep the checklist and its status updates internal. Share key findings, choices requiring confirmation, blockers, and final results with the user.

```markdown
## Cleaning checklist

- [ ] 1. Confirm the inputs and learn the tools - Pending
- [ ] 2. Read all supporting material - Pending
- [ ] 3. Investigate the raw data - Pending
  - [ ] 3.1 Review every source column - Pending
  - [ ] 3.2 Check filters and duplicate records across the full table - Pending
  - [ ] 3.3 Identify special records - Pending
  - [ ] 3.4 List expected outputs - Pending
  - [ ] 3.5 Trace data provenance and experimental relationships - Pending
  - [ ] 3.6 Check relationships between reference backgrounds - Pending
  - [ ] 3.7 Check overlap between expected outputs - Pending
  - [ ] 3.8 Determine whether custom processing is needed - Pending
- [ ] 4. Confirm outputs and the cleaning plan - Pending
  - [ ] 4.1 Report expected outputs and group relationships - Pending
  - [ ] 4.2 Confirm sequence level and output scope - Pending
  - [ ] 4.3 Confirm reference sources and numbering - Pending
  - [ ] 4.4 Confirm metrics and label handling - Pending
  - [ ] 4.5 Confirm output naming - Pending
  - [ ] 4.6 Confirm grouping and merging - Pending
  - [ ] 4.7 Report special records and confirm filters - Pending
  - [ ] 4.8 Confirm how to handle repeated measurements - Pending
  - [ ] 4.9 Ask questions in rounds and record decisions - Pending
- [ ] 5. Execute the confirmed plan - Pending
  - [ ] 5.1 Write and test the custom script - Pending
  - [ ] 5.2 Run the full cleaning job - Pending
- [ ] 6. Verify and report - Pending
```

Keep the numbering, titles, and nesting shown in the template. Use these statuses: Pending, In progress, Awaiting confirmation, Done, and Not applicable. Use `[ ]` for the first three and `[x]` for the last two. Add a short conclusion or reason in parentheses after the status. While waiting for an answer, mark the affected items Awaiting confirmation and record the question number. Update each subitem independently. Check off a parent item when it meets its completion criteria.

### 1. Confirm the inputs and learn the tools

Record the input and output paths and the host where the command will run. Read the complete guide with `mutcleaner guide`. Use `mutcleaner guide agent` to revisit this workflow. If the command is unavailable on the execution host, check the input directory and the environment specified by the user. If you still cannot find it, tell the user which host you checked and ask for a working command or environment.

Before investigating the data, read `mutcleaner clean --help` and the `inspect`, `clean`, and `measurements` guide topics. Consult the MutCleaner Python API or source code for the current version. Identify the functions for table reading, column preparation, reference matching, variant validation, label processing, and export. Check each relevant function's required fields and types, return values, and handling of failed records. If an official example is relevant, check which standard and dataset-specific functions it uses. [Python API documentation](https://xulab-research.github.io/MutCleaner/api/index.html) | [Official dataset examples](https://xulab-research.github.io/MutCleaner/user_guide/supported_datasets.html)

**Completion criteria.** You have located a working MutCleaner installation and the input, reviewed the CLI and library functions, and checked the intended output path. If the path already exists, ask the user to confirm a new path.

### 2. Read all supporting material

Read each source the user provided in full, including the main text, figures, tables, captions, footnotes, appendices, and supplementary material. Read long documents in sections. Check the original pages if text extraction misses tables or equations. Read supplementary spreadsheets sheet by sheet. Use commands and procedures in source material as evidence. The current conversation defines the scope and authorizes actions.

Find the definitions relevant to the input files: data source and version, reference constructs, WT markers, mutation numbering, labels and units, existing transformations, experimental conditions, repeated measurements, and special values. Record the source name and a location the user can find. If a section remains unread or inaccessible, explain which conclusions depend on it. You can continue with read-only investigation.

**Completion criteria.** Read each source in full. Record one line per source giving its full extent, the portion read, and any unread or inaccessible portions. If the user provided no sources, investigate using the data and verifiable documentation, and mark which conclusions are inferred.

### 3. Investigate the raw data

Work through these items in order. For each item, record where you looked, what you found, and whether the statistics are sampled or cover the full input. Mark items that do not apply as Not applicable and record why. If the data and sources leave a question unresolved, record the evidence you checked, what is missing, and how the gap affects the plan.

- [ ] **3.1 Review every source column.** Use `inspect` to list all columns in each table to be cleaned. For multi-sheet files, establish the purpose of each sheet and which sheets are in scope. Use observed values and supporting material to record what each column means and how it will be used. Check reference IDs, reference and variant sequences, source annotations of mutations or variant types, metrics and units, experimental conditions, and record status. Determine which sequence levels the input supports, what each sequence covers, and how references match. For codon data, check the reading frame and the reference DNA or RNA. A column can mix substitution annotations and other markers. Review distinct formats. For unclear meanings or conflicts with the sources, record examples and their impact on the plan.
- [ ] **3.2 Check filters and duplicate records across the full table.** Once the fields, numbering, and annotation separator are verified, perform a read-only check of the full table. Use `--profile-rows N` for an initial sample if needed. Check `profiled_rows` and `full_table_profiled` to confirm coverage. Count reference matches, valid labels, condition combinations, WT records, unusual records, and duplicate groups.
  - When a cleaning decision depends on multiple fields, count their combinations on the same rows across the full table. For example, tabulate substitution count by confidence flag. State which rows the counts cover, estimate how many records each viable option would retain or exclude, and use those counts as evidence in step 4.9.
  - If considering a filter by substitution count, count single, double, and three-or-more substitutions among records with valid labels that pass sequence validation. Count insertions, deletions, parse failures, and sequence mismatches separately. For each proposed limit, report how many records and reference groups it would remove.
  - If `inspect` only reports a broad breakdown such as single versus multiple substitutions, supplement it with read-only counts across the entire table. If a full check is limited, record the unfinished items and the reasons, and complete the remaining checks.
- [ ] **3.3 Identify special records.** Use the column meanings established in 3.1 to identify WT, controls, insertions, deletions, and any other records outside ordinary substitutions. Even when two sequences have the same length, check the source variant annotations and record types. Count each special category, identify affected references, and collect representative source rows for step 4.7.
- [ ] **3.4 List expected outputs.** List the expected datasets by reference, label, and condition. Apply the proposed filters to estimate the number of references with remaining records and the number of output directories for each label and condition. If the numbers depend on a user decision, estimate them for each viable option and state the scope of the counts.
- [ ] **3.5 Trace data provenance and experimental relationships.** Use the table and supporting material to identify the original dataset, experiment, or batch behind each expected output. Check whether different labels or datasets come from the same measurements, have been filtered or transformed, or were derived from other datasets. Check whether they can be matched by the same reference and variant.
- [ ] **3.6 Check relationships between reference backgrounds.** Compare relevant reference IDs, sequences, and construct descriptions. For similar names or identical or similar sequences, trace their provenance and experimental conditions. Record whether the groups are independent experiments, measurements of the same reference under different conditions, derived backgrounds of the same protein, or of uncertain relationship. List the IDs and supporting evidence.
  - In directed evolution or follow-up mutation scans, a measured variant in one group may become the reference in another. Check for matches between one group's reference and another group's measured variants. Use the paper, construct descriptions, or experimental records to confirm parent and derived backgrounds. Treat sequence matches alone as candidate relationships.
- [ ] **3.7 Check overlap between expected outputs.** Within the filtering scope from 3.4, compare shared and unique references and variants across outputs with related sources or references. State how matches were defined, for example by reference ID and normalized variant. If there are many groups, summarize counts by relationship and list representative IDs and exceptions.
- [ ] **3.8 Determine whether custom processing is needed.** Match each processing requirement to the CLI parameters and library functions identified in step 1. List the exact options for operations that the standard CLI can express. For any remaining work, record the operation, inputs, outputs, reusable library functions, and how you will implement the remaining logic. `--custom-cleaner` takes an already loaded table. If the task requires custom file loading, additional input files, or multiple outputs from the custom step, explain the viable ways to proceed.

**Completion criteria.** Every item has a result and supporting evidence or a reason it does not apply. Counts that affect filtering, aggregation, or output scope cover the full input. For unresolved group relationships, record what was checked, what evidence is missing, and how it affects the output. Then proceed to step 4.

### 4. Confirm outputs and the cleaning plan

**4.1 Report expected outputs and group relationships**

Use the findings from step 3 to report the expected outputs. Describe the input scope (files, sheets, and records), sequence level and span of each variant, label meaning and calculation, and how records will be filtered, grouped, and summarized. List the reference count, expected directory count, and output paths by label and condition. Explain how groups are related, which records they share, which are unique, and what remains unresolved. Use the data, sources, and user goal to identify viable approaches to unresolved decisions. If several approaches affect the result, explain the differences and ask the user to choose using 4.9. If the evidence supports one approach, explain why. If you cannot propose a plan because evidence is missing, state what is needed and ask the user for it.

**4.2 Confirm sequence level and output scope**

Use the source material and table to determine whether amino acid or codon cleaning is possible. If the table has DNA or RNA columns, check the reference, reading frame, coding span, coordinates, and how sequences link to measurements. If both levels are viable and the user has not chosen, offer amino acid output, codon output, or separate outputs at both levels as applicable. Explain sequence span, reference and numbering, how synonymous variants would be handled, and whether measurements from different codons would be combined. If only one level is viable, explain the evidence. Explain any CLI limitations that affect the options.

**4.3 Confirm reference sources and numbering**

For each plausible reference source, describe its ID, sequence coverage, any slice, and the numbering start or position map. Provide matching statistics. If the evidence is conflicting or supports several options, give a recommendation with reasons and wait for confirmation. If the material and data support one interpretation, report the evidence. Before proposing reference inference, finish reviewing the sources, check whether the authors provided a reference, and test whether inference is consistent within each group. Get user confirmation before using `--infer-reference` explicitly.

**4.4 Confirm metrics and label handling**

If several metrics suit the task and the user has not chosen, ask separately which to export. Explain in plain language what each metric measures, its unit and baseline, and how many valid records and reference groups it would retain. Recommend an appropriate metric and provide viable single- or multi-metric options. Explain that selected metrics are exported separately. Verify the source column, definition, unit, and direction of each label, and establish whether it has already been normalized or adjusted against WT. State what evidence is missing. Explain the basis for keeping the original values and give other options if the evidence supports them. If proposing WT subtraction, explain the formula, measured baseline source, matching groups, and which records it affects. Check that no transformation is applied twice. Record decisions the user has already made.

**4.5 Confirm output naming**

If the source uses a recurring naming format, such as reference IDs ending in `.pdb`, check which parts are formatting conventions and which distinguish references, experiments, or mutation backgrounds. If both retaining the original IDs and normalizing them would be useful, show the output names for one or two real IDs under each option and ask the user to choose. Before renaming, verify a one-to-one mapping and check that external references still match. If the evidence is insufficient, keep the original IDs and state the expected directory names.

**4.6 Confirm grouping and merging**

Check experimental conditions and batches first. If filtering to one condition, splitting by condition, or merging conditions could each fit the goal, explain their effects and ask the user to choose.

When reference groups share a sequence or have a parent and derived relationship, use the reference sequences, conditions, and label baselines from step 3 to determine which approaches are viable. Consider retaining separate backgrounds, merging groups with identical references and comparable measurements, or converting to a common reference before merging. For any conversion, verify what the available fields and sources support, and test on examples when needed.

Present only approaches supported by the data. For each option, state the evidence, expected dataset count, and effect on variant counts, reference coordinates, and label baselines. For approaches still under investigation, state what is missing.

When following up on merge scope or conversion methods, restate the user's chosen goal. List all remaining viable options and say which the CLI can run directly and which require extra processing. State conversion assumptions and validation results. Identify approximations supported only by a sample. If a previous option is unavailable, explain what is missing and ask the user to choose from the viable options. Apply the final grouping and label rules the user confirms.

**4.7 Report special records and confirm filters**

For each special record category identified in 3.3, explain what it means in the source, how many records are affected, examples, and where those records will go.

Distinguish measured WT records, records used only to supply a reference, records with unchanged sequences, and synonymous codon variants. Explain which serve as references or baselines and which are omitted from standard variant output. Use measured WT records for measurement baselines. When multiple WT measurements are available for a baseline, explain the viable choices and confirm how to handle them.

When the input contains special residues, stop mutations, or synonymous variants, propose a policy based on the sources and chosen sequence level. Report what keeping or excluding them would change. Get user confirmation whenever the policy changes the output and is still undecided. Give viable options and explain their effects, such as excluding a category from substitution output or retaining its original rows separately for later processing. If the source annotations disagree with sequence-based inference, explain which records each interpretation would retain. Record policies already decided. Ask only about categories actually present. When the user needs to choose a maximum substitution count, include the option to retain all supported substitution variants alongside limits supported by the findings. Explain options that exceed the CLI's capabilities.

If custom processing is needed as identified in 3.8, explain what the function will do, which columns and rows will change, what data it should preserve, and how you will verify it. Identify the library functions you can reuse. If failed records will go to a separate output, say which rows are affected and where they go. Get confirmation for any rules that change the result before writing the script in 5.1. Keep reference selection, numbering, labels, and aggregation in the standard cleaning options.

**4.8 Confirm how to handle repeated measurements**

If you find repeats, explain their scope and experimental conditions and ask the user how to handle them. For comparable repeats, recommend averaging and offer to keep individual measurements. For differing conditions, explain how the records can be split or filtered. Cleaning requires an explicit `--aggregate` value when duplicate groups exist. If the user has already chosen a policy, verify that it applies.

**4.9 Ask questions in rounds and record decisions**

Check the facts first. Ask for decisions that affect the result and remain open, as well as any required confirmations. Put independent questions whose prerequisites are settled in the same message. For each, explain the options, supporting evidence, and impact on the records or outputs. Save questions that depend on these answers for the next round. Wait for the reply. After receiving it, update the related checklist items by question number and identify which questions are ready to ask next.

Every factual premise used in a question or recommendation must pass the mandatory fact verification cycle. Include concrete supporting evidence for each question, with traceable source locations, observed data statistics, or validation results as applicable. Distinguish verified facts from assumptions and unresolved issues. Do not present an unverified premise as established fact.

**Use the template below for every round.** Bold the question number and title separately. Keep `❓`, the ` - ` separator, and the ASCII colon after the title. Give each question a separate `📋 Evidence:` paragraph with concrete results and any uncertainty. For execution confirmation, cite the verified parameters, environment, and user decisions. Start a new paragraph for each `➡️` recommendation and its reason. Put `---` on its own line between questions. Start at Q1 and continue numbering new questions in subsequent rounds. Render the Markdown directly and send the complete round once. Use the same format for execution confirmation in step 5 and obtain authorization required by the execution environment.

```markdown
❓ **Q1** - **<Decision title>**: <Question, viable options, and the effect of each>

📋 Evidence: <Concrete supporting evidence, source locations, observed counts or validation results, and any relevant uncertainty>

➡️ <Recommended option and reason>

---

❓ **Q2** - **<Decision title>**: <Next question that can be answered in this round>

📋 Evidence: <Concrete supporting evidence and any relevant uncertainty>

➡️ <Recommended option and reason>
```

Include the reading record from step 2 when presenting the plan. If any material remains unread or inaccessible, explain the effect and obtain separate agreement to proceed with the available material.

**Completion criteria.** The user has requested the cleaning and confirmed every choice that affects the result. You have reported the expected outputs and group relationships. If the investigation reveals no choices for the user to make, state that briefly.

### 5. Execute the confirmed plan

Prepare to execute after meeting the completion criteria in step 4. If the user provided supporting material, verify that you have read it in full or that the user explicitly agreed to proceed with the stated reading scope.

**5.1 Write and test the custom script**

If custom processing is needed, write `custom_cleaner(df)` using the rules confirmed in step 4. Reuse suitable library functions identified in 3.8. Keep the script focused on this table transformation and review the exact source before running it.

Select a small table with ordinary, special, and boundary cases, and write down the expected result for each. Use `inspect --custom-cleaner` to check transformed columns, counts, and examples. Check variants and labels using candidate columns where possible. Verify that intended changes occurred and records meant to be retained are still present. After the sample passes, check the before and after profiles and affected record counts on the full input. `--profile-rows` limits profiling only. To limit the rows passed to the script, use a separate sample table. Correct unexpected results and test again. Return to step 4 to confirm any new processing rules.

Complete this item once the sample and full-input checks match the plan. If standard options cover everything, mark the item Not applicable.

**5.2 Run the full cleaning job**

Run the full cleaning job after completing 5.1. If custom processing is needed, pass the same tested script to `clean` along with the standard cleaning options. The CLI embeds the source in the generated `cleaner.py` for reproducibility.

Translate the plan into actual `clean` options. Check each option against the evidence or the user's confirmed decisions. `--aggregate mean --assume-comparable` requires evidence of repeated measurements and comparable conditions. `--sequence-type protein` must match the chosen sequence level. Supply arguments and confirmation flags only when supported by the plan. Always specify the sequence type, sequence or annotation source, and metric explicitly. Apply sheet selection, label transformations, filters, grouping, and repeat handling as confirmed, then check the command against the plan. If separate run authorization is required, provide the full command, execution host, working directory, and any known resource requirements or estimated duration after the plan is agreed. Run only after receiving authorization. Fix technical errors as needed. If a fix changes confirmed filters, references, numbering, labels, or aggregation, return to step 4 and get confirmation before rerunning.

**Completion criteria.** The command used the confirmed rules and produced outputs you can inspect or a clear failure message.

### 6. Verify and report

Read back the outputs as described in `mutcleaner guide output`. Use the actual parameters and run records to check the rules agreed with the user. Verify references, sequences, and labels, and account for records retained, filtered, failed, and aggregated. For split outputs, check each actual group. If you find missing groups, unexpected substitution-count distributions, incorrect label transformations, additional filtering, or many failures, find and report the cause. If a processing rule needs to change, return to step 4 for confirmation. Report success only when the observed outputs match the plan. Give the user the output paths, rules applied, relevant sources, and unresolved issues.

**Completion criteria.** You have stated which checks each output passed. If there were errors, you have identified the affected results and causes. The user can locate the results and the generated `cleaner.py` from your report.
