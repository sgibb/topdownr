# topdownr 1.1

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

# MSnbase 0.9

## Changes in version 0.99.0
- First public release.
