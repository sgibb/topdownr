xsd <- xml2::read_xml("https://raw.githubusercontent.com/thermofisherlsms/meth-modifications/release/beta-3.1/xsds/Calcium/3.1/MethodModifications.xsd")
devtools::load_all()
.validMsSettings <- .xmlValidMsSettings(xsd)
devtools::use_data(.validMsSettings, internal=TRUE)
