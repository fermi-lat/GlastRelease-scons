<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                version="1.0">
  <xsl:output method="html" encoding="iso-8859-1" />
  <xsl:template match="primary">
    <html>
      <head>
        <title>Primary Constants</title>
      </head>
      <body>
	<xsl:apply-templates select="version" />
        <xsl:apply-templates select="primCategory" />
      </body>
    </html>
  </xsl:template>
  <!-- Do nothing for remainder -->
  
<xsl:template match="version">
  <h2 align="center">
 GLAST Primary Constants, Version 
<xsl:value-of select="@major" />.<xsl:value-of select="@minor" /></h2>
</xsl:template>

<xsl:template match="primCategory">
<h2>Category <font color="#f08000"><xsl:value-of select="@name" />
</font> </h2>
   <xsl:call-template name="overview" />
     <table border="1" cellpadding="3">
     <tr align="left" bgcolor="#c0ffc0">
         <th>Name</th><th>Value</th><th>Description</th>
     </tr>
      <xsl:apply-templates select="prim" />
     </table>
<br /><xsl:text> &#xa0;</xsl:text><br />

</xsl:template>

<xsl:template name="overview">
  <strong><font color="#208020">Overview:</font></strong>
<xsl:text> &#xa0;</xsl:text>
  <xsl:value-of select="overview" disable-output-escaping="yes" />
<br /><xsl:text> &#xa0;</xsl:text><br />
</xsl:template>

<xsl:template match="prim">
  <tr>
  <td><b><xsl:value-of select="@name" /></b></td>

  <td><xsl:value-of select="@value" />
    <xsl:if test = " @uType = 'length' ">
      <xsl:text> &#xa0;</xsl:text>
      <xsl:value-of select="@unitLength" />
    </xsl:if>
</td>
  <td><xsl:value-of select="." /></td>
  </tr>
</xsl:template>
<xsl:template match="derCategory/overview" ></xsl:template>
</xsl:stylesheet>








