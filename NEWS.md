# topdownr 1.29

# topdownr 1.29.1

- Adapt `.calculateFragments` to changes in `PSMatch::calculateFragments`
  introduced in
  [PSMatch:PR19](https://github.com/rformassspectrometry/PSMatch/pull/19).
- Replace partial matching arguments with complete argument names.

# topdownr 1.25

## Changes in version 1.25.1
- New version for Bioc 3.19 (devel)
- Adapt to fragment functions that are moved from `MSnbase` to `PSMatch`.

## Changes in version 1.19.2

# topdownr 1.19

## Changes in version 1.19.2
- Avoid coercion methods that will be deprecated in Matrix 1.4-2.
- Replace `stopifnot` calls with `if (...) stop` calls.

## Changes in version 1.19.1
- New version for Bioc 3.17 (devel)
- Rerun latest `roxygen2`.
- Remove `if (getGeneric(...))` statements in `AllGenerics.R` do avoid problems
  while reloading the package.

# topdownr 1.18
- New version for Bioc 3.14 (release)

# topdownr 1.17

## Changes in version 1.17.3

- Depending on `mzR` 2.27.5.
- Fix unit test for `.readSpectrum` to adapt to new `mzR` 2.27.5.

## Changes in version 1.17.2

- Fix `roxygen2` warnings.

## Changes in version 1.17.1

- New version for Bioc 3.15 (devel)
- Adapt to new `DFrame`.

# topdownr 1.16
- New version for Bioc 3.14 (release)

# topdownr 1.15

## Changes in version 1.15.1
- Add `rmarkdown` to `Suggests:`;
  see https://github.com/yihui/knitr/issues/1864 for details [2021-07-27].

# topdownr 1.13
- New version for Bioc 3.13 (devel)

## Changes in version 1.13.1
- `as(..., "NCBSet")` now treats neutral losses and modifications as bonds as
  well.
- `readTopDownFiles` gains a new argument `customModifications` to allow
  user-defined modifications. Suggestion and first implementation by
  Maša Babović <masab@bmb.sdu.dk> [2021-03-15].

# topdownr 1.12
- New version for Bioc 3.12 (release)

# topdownr 1.11
- New version for Bioc 3.12 (devel)

# topdownr 1.10
- New version for Bioc 3.11 (release)

# topdownr 1.9
- New version for Bioc 3.11 (devel)

## Changes in version 1.9.4
- Fix unit tests that check for "matrix" class.
  (`class(m)` now returns `c("matrix", "array")` in r-devel) [2019-12-17].

## Changes in version 1.9.3
- Adapt analysis vignette to changed `condition` argument [2019-11-22].

## Changes in version 1.9.2
- Set `conditions="ScanDescription"` as new default for `readTopDownFiles`.
  The creation of *FilterString* IDs in the method files was deprecated since
  over a year. `conditions="FilterString"` is still possible for
  backward-compatibility [2019-11-22].
- Fix error message handling in `.validFilename` and
  `.translateThermoIdToScanId` [2019-11-22].

## Changes in version 1.9.1
- Remove `defaultMs1Settings` and `defaultMs2Settings` [2019-11-18].
- Fix `.rbind` for lists with mixed `data.frame` and `DataFrame` [2019-11-18].

# topdownr 1.8
- New version for Bioc 3.10 (release)

# topdownr 1.7
- New version for Bioc 3.10 (devel)

## Changes in version 1.7.2
- Add `readTopDownSet(..., conditions="ScanDescription")` as a new way to
  read scan conditions (see [#80](https://github.com/sgibb/topdownr/issues/80)/[#81](https://github.com/sgibb/topdownr/issues/81)) [2019-08-08].

## Changes in version 1.7.1
- Remove NEWS file (just keep NEWS.md).
- Never remove `"AgcTarget"` column from `colData` `DataFrame`.
- Strip white spaces from ScanHeadsman output.
- Defunct `defaultMs1Settings` and `defaultMs2Settings`. They will be
  removed in 3.11 [2019-06-19].

# topdownr 1.6
- New version for Bioc 3.9 (release)

# topdownr 1.5
- New version for Bioc 3.9 (devel)

## Changes in version 1.5.6
- Revert changes for NULL indices of `DataFrame` introduced in 1.5.4
  (a419f59, c4bfc1c) because they are fixed upstream in `S4Vectors`.
  Keep unit tests in place. [2019-03-27]

## Changes in version 1.5.5
- Depends on R >= 3.5.0 now, because the seralization format changed in R.

## Changes in version 1.5.4
- Fix for internal `.makeRowNames`/`.colsToLogical`/`.colsToRle`
  on `DataFrame` without any `numeric`/`character` columns.
- Fix unit test that uses `set.seed` (order changed during R-devel upgrade).

## Changes in version 1.5.3
- biocViews: ImmunoOncology added by
  Kayla-Morrell <kayla.morrell@roswellpark.org> [2019-01-04].

## Changes in version 1.5.2
- Add `expandMs1Conditions`, `expandTms2Conditions`,
  `createExperimentsFragmentOptimisation` functions to allow more flexibility
  in method creation; see also [#77](https://github.com/sgibb/topdownr/issues/77) [2018-12-07].
- Modify interface/arguments of `writeMethodXmls` to adapt to new method
  creation workflow (the old interface will be defunct in Bioconductor 3.10 and
  removed in 3.11) [2018-12-07].
- Adapt the `data-generation` vignette to the new workflow [2018-12-07].
- Deprecated `defaultMs1Settings` and `defaultMs2Settings`. They will be
  defunct in Bioconductor 3.10 and removed in 3.11 [2018-12-07].

## Changes in version 1.5.1
- `readTopDownFiles` gains the argument `conditions` to control wheter
  "FilterStrings" or a given number of conditions is used to create condition
  IDs; see [#77](https://github.com/sgibb/topdownr/issues/77) [2018-11-07].

# topdownr 1.4
- New version for Bioc 3.8 (release)

# topdownr 1.3

## Changes in version 1.3.6
- Add Pavel's and Ole's ORCID to DESCRIPTION [2018-10-23].

## Changes in version 1.3.5
- Fix format of roxygen links to foreign packages to avoid link warning in
  `R CMD check` [2018-10-10].

## Changes in version 1.3.4
- Add inst/CITATION file [2018-09-26].

## Changes in version 1.3.3
- Revert commit c6e8dfd: "Adapt to `MSnbase 2.7.2` with internal fragments; see
    [#82](https://github.com/lgatto/MSnbase/issues/82) [2018-06-03]."

## Changes in version 1.3.2
- Use `BiocManager::install` [2018-07-16].

## Changes in version 1.3.1
- Adapt to `MSnbase 2.7.2` with internal fragments; see
    [#82](https://github.com/lgatto/MSnbase/issues/82) [2018-06-03].
- Fix `FragmentViews` start/end/width and labels for internal fragments
    [2018-06-03].
- Fix `as(tds, "MSnSet")` unit test [2018-07-06].
- Use `elementMetadata(..., use.names=FALSE)` in
  `combine,FragmentViews,FragmentViews-method` to avoid duplicated rownames in
  elementMetadata slot [2018-07-06].

## Changes in version 1.3.0
- New version for Bioc 3.8 (devel)

# topdownr 1.2

## Changes in version 1.2.0
- New version for Bioc 3.7 (release)

# topdownr 1.1

## Changes in version 1.1.7
- Add `mz,FragmentViews-method` [2018-02-01].
- Remove internal `fragmentMass` and `fragmentNames` functions [2018-02-22].
- Parse "spectrumId" column of the mzML header to find the scan number (instead
  of the "acquisitionNum") because ProteomDiscover generates non-standard
  "spectrumId" and proteowizard fails to translated it into a valid
  "acquisitionNum". See [#73](https://github.com/sgibb/topdownr/issues/73) for details [2018-02-22].
- Recalculate TotIonCurrent in the main loop of `.readMzMl` [2018-02-22].
- Add `FragmentCoverage` and `BondCoverage` columns to
  `bestConditions,NCBSet-method` [2018-02-23].
- Use retention times to test for correct matching between ScanHeadsman .txt
  output and mzML files; closes [#74](https://github.com/sgibb/topdownr/issues/74); [2018-02-23].

## Changes in version 1.1.6
- Rotate fragment labels (vertical orientation) in `plot` [2018-01-17].
- Replace signature for `updateMedianInjectionTime,TopDownSet-method` to
  `updateMedianInjectionTime,AbstractTopDownSet-method`; closes [#69](https://github.com/sgibb/topdownr/issues/69); see
  also [#71](https://github.com/sgibb/topdownr/issues/71) [2018-01-27].
- Fix `.matchFragments` for `length(fmass) == 0` [2018-01-27].
- Just plot fragments that are present in current `TopDownSet` see [#70](https://github.com/sgibb/topdownr/issues/70)
  [2018-01-27].
- Add `combine,FragmentViews,FragmentViews-method` [2018-01-27].
- Allow to `combine` `TopDownSet` objects with different fragment types;
  closes [#71](https://github.com/sgibb/topdownr/issues/71) [2018-01-27].
- Add `all.equal` for `AbstractTopDownSet` objects [2018-01-27].
- Allow the user to decide how to handle redundant fragment matching. Current
  default is `redundantFragmentMatch="remove"` and
  `redundantIonMatch="remove"`. This will reduce the number of fragment
  matches. Choose `"closest"` for both to get the old behaviour.
  See also [#72](https://github.com/sgibb/topdownr/issues/72) [2018-01-29].
- `TopDownSet` object store the matching `tolerance` and strategies
  (`redundantIonMatch`, `redundantFragmentMatch`). `AbstractTopDownSet` and
  `NCBSet` lost their `tolerance` slot. Saved objects need to be recreated
  [2018-01-30].
- `bestConditions,NCBSet-method` returns a 5-column matrix now.
  Colums are: Index, FragmentsAddedToCombination, BondsAddedToCombination,
  FragmentsInCondition, BondsInCondition; see [#52](https://github.com/sgibb/topdownr/issues/52) [2018-01-30].

## Changes in version 1.1.5
- Keep full filename (before `basename` was used) in `AbstractTopDownSet`
  objects [2017-12-28].
- Add `plot,TopDownSet-method` [2017-12-29].
- `bestConditions,NCBSet-method` gains a new argument `maximise` that allows to
  optimise for number of fragments or bonds covered (default: `"fragments"`);
  see [#52](https://github.com/sgibb/topdownr/issues/52) [2018-01-15].

## Changes in version 1.1.4
- Add missing export of `combine` and documentation [2017-12-28].
- Resave `tds` example data set to reflect changes in `colData` introduced in
  version 1.1.2 [2017-12-28].

## Changes in version 1.1.3
- Add `conditionNames,AbstractTopDownSet-method` to access
  `rownames(colData(tds))` [2017-12-23].
- Add `updateConditionNames,AbstractTopDownSet-method`
  (closes [#60](https://github.com/sgibb/topdownr/issues/60)) [2017-12-23].
- Turn `updateMedianInjectionTime,TopDownSet-method` into
  `updateMedianInjectionTime,AbstractTopDownSet-method` to work with
  `TopDownSet` and `NCBSet` objects [2017-12-27].
- Add `combine,AbstractTopDownSet-method` to combine multiple
  `TopDownSet`/`NCBSet` objects (closes [#69](https://github.com/sgibb/topdownr/issues/69)) [2017-12-28].

## Changes in version 1.1.2
- Add `.rbind` to combine scan and method information with different number of
  colums (could happen when CID/HCD and UVPD scans are taken independently with
  different software versions) [2017-12-22].
- Don't replace NA values with zeros in the `colData` [2017-12-22].
- Convert On/Off `character` columns in scan and method information to
  `logical` [2017-12-22].
- Fix `.camelCase` to avoid "TIC" to "TIc" and "UseCalibratedUVPDTimeMs2" to
  "UseCalibrateduvpdTimems2" conversion (now: "Tic" and
  "UseCalibratedUvpdTimeMs2") [2017-12-22].

## Changes in version 1.1.1
- Respect assigned intensity in conditions for `bestConditions,NCBSet-method`
  and `fragmentationMap`
  (closes [#62](https://github.com/sgibb/topdownr/issues/62)) [2017-12-02].
- Fix explanation of random forest barchart in analysis vignette [2017-12-02].
- Create all fragmentation methods in `.readScanHeadsTable` to avoid error if
  any is missing (fixes [#68](https://github.com/sgibb/topdownr/issues/68)) [2017-12-20].
- Never remove Activation column in `colData` (even not if
  `readTopDownFiles(..., dropNonInformativeColumns=TRUE)`) [2017-12-20].
- Allow UVPD in `fragmentationMap,NCBSet-method` [2017-12-20].
- Add new method: `updateMedianInjectionTime,TopDownSet-method`
  (closes [#66](https://github.com/sgibb/topdownr/issues/66)) [2017-12-20].

## Changes in version 1.1.0
- New version for Bioc 3.7 (devel)

# topdownr 1.0

## Changes in version 1.0.0
- New version for Bioc 3.6 (release)

# topdownr 0.9

## Changes in version 0.99.0
- First public release.
