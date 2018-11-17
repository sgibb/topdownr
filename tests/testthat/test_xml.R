context("xml")

xsd <- '<?xml version="1.0" encoding="utf-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema" >
  <xs:simpleType name="Version">
    <xs:restriction base="xs:string">
      <xs:enumeration value="1"/>
      <xs:enumeration value="2"/>
    </xs:restriction>
  </xs:simpleType>
  <xs:element name="FullMSScan">
    <xs:complexType>
      <xs:all>
        <xs:element name="FirstMass" type="xs:double" minOccurs="0" maxOccurs="1" />
        <xs:element name="LastMass" type="xs:double" minOccurs="0" maxOccurs="1" />
        <xs:element name="Microscans" type="xs:int" minOccurs="0" maxOccurs="1" />
      </xs:all>
    </xs:complexType>
  </xs:element>
  <xs:element name="TMSnScan">
    <xs:complexType>
      <xs:all>
        <xs:element name="AgcTarget" type="xs:double" minOccurs="0" maxOccurs="1" />
        <xs:element name="ETDReactionTime" type="xs:double" minOccurs="0" maxOccurs="1" />
        <xs:element name="CIDCollisionEnergy" type="xs:double" minOccurs="0" maxOccurs="1"/>
        <xs:element name="CIDActivationQ" type="xs:double" minOccurs="0" maxOccurs="1"/>
        <xs:element name="Ms2CIDCollisionEnergy" type="xs:double" minOccurs="0" maxOccurs="1"/>
        <xs:element name="Ms2CIDActivationQ" type="xs:double" minOccurs="0" maxOccurs="1"/>
        <xs:element ref="MassList" minOccurs="0" maxOccurs="1"/>
        <xs:element name="ScanDescription" minOccurs="0" maxOccurs="1"/>
      </xs:all>
    </xs:complexType>
  </xs:element>
</xs:schema>'

test_that(".xmlChildrenNameType", {
    skip_if_not_installed("xml2")
    expect_equal(
        .xmlChildrenNameType(xml2::read_xml(xsd), "FullMSScan"),
        cbind(
            name=c("FirstMass", "LastMass", "Microscans"),
            class=c("double", "double", "integer")
        )
    )
})

test_that(".xmlEnumeration", {
    skip_if_not_installed("xml2")
    expect_equal(
        .xmlEnumeration(xml2::read_xml(xsd), "Version"),
        c("1", "2")
    )
})

test_that(".xmlTmsnScanParameters", {
    skip_if_not_installed("xml2")
    expect_equal(
        .xmlTmsnScanParameters(xml2::read_xml(xsd)),
        cbind(
            name=c(
                "AgcTarget", "ETDReactionTime", "CIDCollisionEnergy",
                "CIDActivationQ", "MassList", "ScanDescription"
            ),
            class=c(rep("double", 4), "MassListRecord", "character")
        )
    )
})

test_that(".xmlValidMsSettings", {
    skip_if_not_installed("xml2")
    expect_equal(
        .xmlValidMsSettings(xml2::read_xml(xsd)),
        cbind(
            name=c(
                "FirstMass", "LastMass", "Microscans", "AgcTarget",
                "ETDReactionTime", "CIDCollisionEnergy", "CIDActivationQ",
                "MassList", "ScanDescription"
            ),
            class=rep(
                c("double", "integer", "double", "MassListRecord", "character"),
                c(2, 1, 4, 1, 1)
            ),
            type=rep(c("MS1", "TMS2", "ETD", "CID", "TMS2"), c(3, 1, 1, 2, 2))
        )
    )
})
