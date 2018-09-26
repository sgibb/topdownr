xsd <- xml2::read_xml("https://raw.githubusercontent.com/thermofisherlsms/meth-modifications/release/beta-3.1/xsds/Calcium/3.1/MethodModifications.xsd")
devtools::load_all()

.validMsSettingsCalcium3.1 <- .xmlValidMsSettings(xsd)

devtools::use_data(
    .validMsSettingsCalcium3.1,
    internal=TRUE,
    overwrite=TRUE
)
