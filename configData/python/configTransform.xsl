<?xml version="1.0"?>
<!DOCTYPE xsl:stylesheet [
<!ENTITY nbsp "&#160;">
]>

<xsl:stylesheet 
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:f="http://xmlns.4suite.org/ext"
    version="1.0">

  <xsl:output method="html"/>

  <xsl:template match="/">
      <xsl:apply-templates select="ConfigReport"/>
      <xsl:apply-templates select="PrecinctReport"/>
      <xsl:apply-templates select="HeritageReport"/>
  </xsl:template>

  <xsl:template match="ConfigReport">
    <HTML>
      <HEAD>
        <TITLE>Configuration report for MOOT config name 
              <xsl:value-of select="@ConfigName"/> (Key =  
              <xsl:value-of select="@ConfigKey"/>)
        </TITLE>
        <STYLE media ="print" type ="text/css">#toc {display:none}p {font: 12pt serif; page-break-inside:avoid}h1,h2 {page-break-before:always}</STYLE>
        <STYLE media ="screen" type ="text/css">h1,h2 {border-top-style:solid}</STYLE>
      </HEAD>
      <BODY>
        <h1>Configuration report for MOOT config name 
              <xsl:value-of select="@ConfigName"/> (Key =  
              <xsl:value-of select="@ConfigKey"/>)
        </h1>
        <p>&nbsp;</p>
        <p>Configuration key <xsl:value-of select="@ConfigKey"/> was compared
           to the baseline configuration <xsl:value-of select="@BaselineKey"/>
           by user <xsl:value-of select="@User"/> on  <xsl:value-of select="@Date"/>
        </p>
        <hr/>
        <xsl:apply-templates select="*"/>
      </BODY>
    </HTML>
  </xsl:template>


  <xsl:template match="PrecinctReport">
    <HTML>
      <HEAD>
        <TITLE>Precinct report for <xsl:value-of select="@Name"/> 
               (Vote key=<xsl:value-of select="@VoteKey"/>, 
               Alias=<xsl:value-of select="@Alias"/>)
        </TITLE>
        <STYLE media ="print" type ="text/css">#toc {display:none}p {font: 12pt serif; page-break-inside:avoid}h1,h2 {page-break-before:always}</STYLE>
        <STYLE media ="screen" type ="text/css">h1,h2 {border-top-style:solid}</STYLE>
      </HEAD>
      <BODY>
        <h1>Precinct report for <xsl:value-of select="@Name"/>
               (Vote key=<xsl:value-of select="@VoteKey"/>, 
               Alias=<xsl:value-of select="@Alias"/>)
        </h1>
        <p>&nbsp;</p>
        <p>Precinct <xsl:value-of select="@Name"/> from MOOT configuration key
           <xsl:value-of select="@ConfigKey"/> was compared
           to baseline configuration <xsl:value-of select="@BaselineKey"/>
        </p>
        <hr/>
        <xsl:apply-templates/>
      </BODY>
    </HTML>
  </xsl:template>

  <xsl:template match="HeritageReport">
    <HTML>
      <HEAD>
        <TITLE>Heritage report for Configuration <xsl:value-of select="@ConfigKey"/> 
        </TITLE>
        <STYLE media ="print" type ="text/css">#toc {display:none}p {font: 12pt serif; page-break-inside:avoid}h1,h2 {page-break-before:always}</STYLE>
        <STYLE media ="screen" type ="text/css">h1,h2 {border-top-style:solid}</STYLE>
      </HEAD>
      <BODY>
        <h1>Heritage report for Configuration <xsl:value-of select="@ConfigKey"/> 
        </h1>
        <hr/>
        <xsl:apply-templates/>
      </BODY>
    </HTML>
  </xsl:template>


  <xsl:template match="ConfigInfo">
    <h3>Configuration database information</h3>
    <table border='1'>
      <tr><td align='left'>Name/Alias:</td>
          <td align='left'><xsl:value-of select="@ConfigName"/></td></tr>
      <tr><td align='left'>DB Key:</td>
          <td align='left'><xsl:value-of select="@Key"/></td></tr>
      <tr><td align='left'>Creation Date:</td>
          <td align='left'><xsl:value-of select="@Date"/></td></tr>
      <tr><td align='left'>Active:</td>
          <td align='left'><xsl:value-of select="@Active"/></td></tr>
      <tr><td align='left'>Status:</td>
          <td align='left'><xsl:value-of select="@Status"/></td></tr>
      <tr><td align='left'>Description:</td>
          <td align='left'><xsl:value-of select="@Description"/></td></tr>
      <xsl:if test="@VoteKey">
        <tr><td align='left'>Vote Key:</td>
            <td align='left'><xsl:value-of select="@VoteKey"/></td></tr>
      </xsl:if>
      <xsl:if test="@FileName">
        <tr><td align='left'>Vote File:</td>
            <td align='left'><xsl:value-of select="@FileName"/></td></tr>
      </xsl:if>
    </table>
  </xsl:template>

  <xsl:template match="BaselineInfo">
    <h3>Baseline configuration database information</h3>
    <table border='1'>
      <tr><td align='left'>Name/Alias:</td>
          <td align='left'><xsl:value-of select="@ConfigName"/></td></tr>
      <tr><td align='left'>DB Key:</td>
          <td align='left'><xsl:value-of select="@Key"/></td></tr>
      <tr><td align='left'>Creation Date:</td>
          <td align='left'><xsl:value-of select="@Date"/></td></tr>
      <tr><td align='left'>Active:</td>
          <td align='left'><xsl:value-of select="@Active"/></td></tr>
      <tr><td align='left'>Status:</td>
          <td align='left'><xsl:value-of select="@Status"/></td></tr>
      <tr><td align='left'>Description:</td>
          <td align='left'><xsl:value-of select="@Description"/></td></tr>
      <xsl:if test="@VoteKey">
        <tr><td align='left'>Vote Key:</td>
            <td align='left'><xsl:value-of select="@VoteKey"/></td></tr>
      </xsl:if>
      <xsl:if test="@FileName">
        <tr><td align='left'>Vote File:</td>
            <td align='left'><xsl:value-of select="@FileName"/></td></tr>
      </xsl:if>
    </table>
  </xsl:template>

  <xsl:template match="PrecinctInfo">
    <h3>Vote database information</h3>
    <table border='1'>
      <tr><td align='left'>Precinct:</td>
          <td align='left'><xsl:value-of select="@PrecinctName"/></td></tr>
      <tr><td align='left'>Alias:</td>
          <td align='left'><xsl:value-of select="@Alias"/></td></tr>
      <tr><td align='left'>Vote Key:</td>
          <td align='left'><xsl:value-of select="@Key"/></td></tr>
      <tr><td align='left'>Creation Date:</td>
          <td align='left'><xsl:value-of select="@Date"/></td></tr>
      <tr><td align='left'>Creator:</td>
          <td align='left'><xsl:value-of select="@User"/></td></tr>
      <tr><td align='left'>Status:</td>
          <td align='left'><xsl:value-of select="@Status"/></td></tr>
      <tr><td align='left'>Description:</td>
          <td align='left'><xsl:value-of select="@Description"/></td></tr>
      <tr><td align='left'>Vote File:</td>
          <td align='left'><xsl:value-of select="@FileName"/></td></tr>
    </table>
    <xsl:if test="child::AncillaryInfo">
      <h4>Ancillary files used:</h4>
      <xsl:for-each select="AncillaryInfo">
        Type:<xsl:value-of select="@Name"/>,
        Path:<xsl:value-of select="@FileName"/>,
        Key:<xsl:value-of select="@Key"/>,
        Created by:<xsl:value-of select="@User"/> on <xsl:value-of select="@Date"/>,
        Status:<xsl:value-of select="@Status"/>
        <br/>
      </xsl:for-each>
    </xsl:if>
    <xsl:if test="child::ParamInfo">
      <h4>Parameter files produced:</h4>
      <xsl:for-each select="ParamInfo">
        Type:<xsl:value-of select="@Name"/>,
        Path:<xsl:value-of select="@FileName"/>,
        Key:<xsl:value-of select="@Key"/>,
        Created by:<xsl:value-of select="@User"/> on <xsl:value-of select="@Date"/>,
        Status:<xsl:value-of select="@Status"/>
        <br/>
      </xsl:for-each>
    </xsl:if>

  </xsl:template>

  <xsl:template match="Precincts">
    <hr/>
    <h4>Precinct reports:</h4>
    <xsl:apply-templates select="*"/>
  </xsl:template>

  <xsl:template match="PrecinctReportLink">
    <a>
      <xsl:attribute name='href'><xsl:value-of select="@FileLink"/>.html</xsl:attribute>
      <xsl:value-of select="@Name"/>
    </a>
    <br/>
  </xsl:template>

  <xsl:template match="ReportSection">
    <hr/>
    <h4><xsl:value-of select="@Name"/></h4>
    <xsl:apply-templates/>
    <hr/>
  </xsl:template>

  <xsl:template match="Intent">
    <h5>Intent:</h5>
    <table><tr><td><xsl:value-of select="."/></td></tr></table>
  </xsl:template>

  <xsl:template match="TextInclude">
    <h6>Included File:&nbsp;&nbsp;<xsl:value-of select="@FileName"/></h6>
    <xsl:choose>
      <xsl:when test="@NoConvert='True'">
        <pre>
          <xsl:apply-templates select="TextLine"/>
        </pre>
      </xsl:when>
      <xsl:otherwise>
        <p>
          <xsl:apply-templates select="TextLine"/>
        </p>
      </xsl:otherwise>
    </xsl:choose>
    <xsl:if test="@NLines">
      <em><b>Note:</b>File truncated to <xsl:value-of select="@NLines"/> lines</em><br/>      
    </xsl:if>
    <p><a>
      <xsl:attribute name='href'><xsl:value-of select="@FileName"/></xsl:attribute>
      Full Output
    </a></p>
  </xsl:template>

  <xsl:template match="TextLine">
    <xsl:choose>
      <xsl:when test="@NoConvert='True'">
        <xsl:value-of select="@Text"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="@Text" disable-output-escaping='yes'/>
      </xsl:otherwise>
    </xsl:choose>

  </xsl:template>

  <xsl:template match="Image">
    <p align="center"><b><xsl:value-of select="@Title"/></b></p>
    <p align="center"><img>
      <xsl:attribute name="src"><xsl:value-of select="@FileName"/></xsl:attribute>
      <xsl:attribute name="alt"><xsl:value-of select="@FileName"/></xsl:attribute>
      <xsl:if test="@Height">
        <xsl:attribute name="height"><xsl:value-of select="@Height"/></xsl:attribute>
      </xsl:if>
      <xsl:if test="@Width">
        <xsl:attribute name="width"><xsl:value-of select="@Width"/></xsl:attribute>
      </xsl:if>
    </img></p>
    <p align="center"><em><xsl:value-of select="."/></em></p>
  </xsl:template>

  <xsl:template match="Comment">
    <h5>Comment by <xsl:value-of select="@User"/> at <xsl:value-of select="@Date"/>:</h5>
    <blockquote><xsl:value-of select="."/></blockquote>
  </xsl:template>

  <xsl:template match="LinkTo">
    <p><a>
      <xsl:attribute name="href">
        <xsl:choose>
          <xsl:when test="@NoConvert='False'">
            <xsl:value-of select="f:replace('.xml', '.html', @FileName)"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="@FileName"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:attribute>
      <xsl:value-of select="@Description"/>
    </a></p>
  </xsl:template>

  <xsl:template match="*">
    <xsl:element name="{name(.)}">
      <xsl:for-each select="@*">
        <xsl:attribute name="{name(.)}">
            <xsl:value-of select="."/>
        </xsl:attribute>
      </xsl:for-each>
      <xsl:apply-templates/>
    </xsl:element>
  </xsl:template>


</xsl:stylesheet>
