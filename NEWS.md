# topdownr 1.1

## Changes in version 1.1.2
- Add `.rbind` to combine scan and method information with different number of
  colums (could happen when CID/HCD and UVPD scans are taken independently with
  different software versions) [2017-12-22].
- Don't replace NA values with zeros in the `colData` [2017-12-22].

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
