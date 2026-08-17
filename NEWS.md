# kinfitr 0.9.3

## BIDS parsing

* **New parser: `bids_parse_filenames()` enumerates the acquisitions on disk.**
  The previous parser worked by completion: it expanded the combinations of
  every entity value observed anywhere in the study, merged each file into the
  combinations it was compatible with, and dropped combinations with nothing to
  merge — substituting defaults (`ses-01`, `task-rest`, `run-01`, and
  placeholder values for `trc`, `rec` and `acq`) for entities a filename did
  not carry. On simple, uniform datasets this resolved itself correctly. On
  complex ones — subjects with and without sessions in one study, multiple runs
  or reconstructions, entities used by some scans and not others — it could
  accidentally create acquisitions that do not exist, and assign a scan entity
  values its own files never carried. That second effect made a sessionless
  subject take on its neighbours' `ses-test`, after which its own files failed
  to join against that identity and the subject silently dropped out of results
  whenever it was analysed alongside others.

  The new parser returns one row per acquisition actually on disk — any
  `_pet.nii`, `_pet.nii.gz` or `_pet.json` inside a `pet/` directory, the image
  and its sidecar counting once — carrying exactly the entities its filename
  and directory path name. Nothing is substituted: a study that does not use
  `trc` has no `trc` column.

* **Files attach to acquisitions by the subset rule.** An entity a file names
  must match the acquisition; an entity it omits imposes no constraint. Blood
  naming no `rec` therefore reaches both of a subject's reconstructions, while
  blood claiming `rec-A` never attaches to anything else.

* **`bids_parse_files()` is deprecated.** It keeps its old behaviour unchanged,
  now warns on use, and will be removed in 2027. `bids_parse_study()` keeps its
  name and contract — parse a study, get its measurements — built on the new
  enumeration.

## Sidecar clashes no longer cost the whole study

* **A dataset that breaks the inheritance principle now costs one acquisition,
  not the parse.** The BIDS spec is explicit that there MUST NOT be more than
  one metadata file applicable to a data file at a single directory level, and
  `bids_resolve_sidecars()` refuses such a dataset rather than picking by
  enumeration order. That error previously propagated out of
  `bids_parse_study()`, so one ambiguous sidecar took every unrelated subject
  in the study down with it. It is now reported per acquisition — naming both
  sidecars and the measurement excluded — and the rest of the study parses,
  matching how `bids_associate_blood()` already handled the equivalent case.
  Worth knowing: the BIDS validator does not currently flag this layout, so
  kinfitr may be the only thing that tells you about it.

## Plotting

* **`qplot()` replaced throughout**, which is deprecated in ggplot2 3.4.0 and
  was warning on every t\* plot.

* **Scale limits that silently dropped data are now coordinate limits.** The
  t\* and fit plots set their view with `xlim()`/`ylim()`, which removes points
  outside the range before drawing — so a `geom_line` through them broke at the
  gap, and every panel warned. They now use `coord_cartesian()`, which clips
  the view and keeps the data.

* **The undefined first point of a linearised plot is declared as such.** Logan
  and Patlak both divide by a concentration that is zero at *t* = 0, so their
  first fitted value is `NA` by construction. Those layers now pass
  `na.rm = TRUE` rather than warning once per panel about a value that is never
  going to be there.

* `get_units_radioactivity()` is **exported**. petfit and bloodstream were
  reaching it through `:::`.


## Sidecar inheritance

* **New `bids_resolve_sidecars()`: inherited metadata is merged
  deterministically.** When several sidecars apply to a data file they are
  merged field by field from the study root inwards, a nearer file overriding a
  further one. Previously the first applicable sidecar in filesystem
  enumeration order won, so inherited fields — frame times included — could
  depend on the order files happened to be listed in. Two applicable sidecars
  in the same directory are now an error rather than an arbitrary pick.

* **Inheritance sidecars are metadata, not acquisitions.** A `_pet.json`
  serving several images — `sub-01_pet.json` beside `run-01` and `run-02` — is
  their shared metadata rather than an additional run-less acquisition. A
  `_pet.json` matching no image still counts as a measurement, so an
  acquisition whose image is not found (most commonly, blood collected before
  the image is available) remains usable.


## Blood association

* **New `bids_associate_blood()`: the blood-to-acquisition mapping is explicit
  and checked.** Blood must sit in the same `pet/` directory as the
  acquisition's own PET data — not merely some `pet/` directory, so
  subject-level blood cannot silently attach to a session-level acquisition —
  and must carry a `recording` entity (both spec requirements), with at most
  one tsv/json pair per recording per acquisition. Several recordings per
  acquisition — manual samples alongside an autosampler — work as before;
  two files claiming the same recording is an error rather than an arbitrary
  pick.

* **A tsv and its json must agree to be a pair.** Every entity the json names
  must be named by the tsv with the same value (naming fewer is fine, per
  inheritance); a `rec-A` json beside a tsv naming no `rec` may describe
  different data, and is warned about rather than paired by their shared
  `recording` label.

* **Unusable blood is reported, not silently absent.** A recording without a
  complete tsv/json pair is named in a warning — including when another
  recording of the same acquisition is complete — instead of being
  indistinguishable from blood that was never collected. Recording labels
  other than `manual` and `autosampler` are likewise named rather than dropped
  without comment.

* Blood no longer leaks across sessions: a session with no blood files gets
  none, rather than borrowing a neighbour's.

## Robustness

* **One malformed acquisition no longer aborts the study parse.** A missing or
  unreadable `_pet.json`, or one lacking frame timing, previously failed the
  whole parse with `object 'dur' not found`; it now excludes that one
  acquisition with a warning naming it and the reason.

* Entity labels differing only in case (`trc-PF974` / `trc-pf974`) produce a
  warning pointing at the BIDS validator — they are treated as distinct, which
  is rarely what was meant.

* `+` in an entity label is accepted. It is legal in BIDS, and previously
  caused the whole key-value pair to be dropped, so `task-A+B` parsed as no
  task at all.

* Two images resolving to the same acquisition is an explicit error rather
  than a silent pick of whichever came first.

## Derivatives

* **New `bids_parse_derivatives()`: derived-data folders get their own
  parser.** A derivatives folder contains no PET images, so the raw-study
  parser was always the wrong tool for it — the completion behaviour fabricated
  entities there too. The new parser returns one row per derived file with
  three keys: `source_key` (which acquisition it derives from, the directory
  completing entities the filename omits), `artifact_key` (which product — a
  tsv and its json sidecar share one), and `analysis_scope_key` (which
  analysis, for group-level files). Group-level files get no `source_key`, so
  they cannot be joined to an acquisition by accident, and the same acquisition
  analysed in two side-by-side analysis folders yields distinct artifacts.

* **The bloodstream importers use it**, fixing three real failures: entity
  columns filled with the old substituted defaults, which then broke joins
  against real TAC data; analysis-level `bloodstream_config.json` files read as
  per-measurement AIF fits; and a crash on an input function missing its
  sidecar, which is now named in a warning and skipped.
