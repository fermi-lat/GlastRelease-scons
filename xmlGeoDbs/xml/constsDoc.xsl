<?xml version="1.0"?>
<!-- $Id$  -->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                version="1.0">
  <xsl:output method="html" encoding="iso-8859-1" />
  <xsl:template match="/">
    <xsl:apply-templates match = "gdd" />
  </xsl:template>
  <xsl:template match="gdd">
    <html>
      <head>
        <title>Detector Constants</title>
      </head>
      <body>
      <h1 align="center">Primary and Derived Constants</h1>
     <xsl:call-template name="catList" />
     <xsl:apply-templates select="constants" />
     <xsl:apply-templates select="section" />
<!-- add our own source information -->
<h3 align="center">Genesis</h3>
     <p>
      Derived from xml source 
      <xsl:value-of 
         select="substring-before(substring-after(@CVSid,'$'),'$')" /> 
<!--      <xsl:value-of 
         select="@CVSid" /> <br />   -->
      created by constsDoc.xsl, <br />
      <xsl:value-of 
         select="substring-before(substring-after('$Id$','$'),'$')" /> 
     </p>
     <xsl:apply-templates select="source" />
      </body>
    </html>
  </xsl:template>

  <xsl:template name="catList">
    <xsl:for-each select="constants/primary">
      <h3>Primary Constants</h3>
      <ul>
      <xsl:for-each select="primCategory">
        <xsl:call-template name="listCat" />
      </xsl:for-each>
      </ul>
    </xsl:for-each>
    <xsl:for-each select="constants/derived">
      <h3>Derived Constants</h3>
      <ul>
      <xsl:for-each select="derCategory">
        <xsl:call-template name="listCat" />
      </xsl:for-each>
      </ul>
    </xsl:for-each>
  </xsl:template>

  <xsl:template name="listCat" >
    <li><a href="#{@name}" > <xsl:value-of select="@name" /> </a> </li>
  </xsl:template>

  <xsl:template match="source">
     <p>..in turn derived from 
     <xsl:value-of select="@inputId" />, <br /> created by
     <xsl:value-of select="@creatorId" />
     </p>
  </xsl:template>

  <xsl:template match="constants" >
     <xsl:apply-templates select="primary" />
     <xsl:apply-templates select="derived" />
  </xsl:template>

  <xsl:template match="primary">
	<xsl:apply-templates select="version" />
        <xsl:apply-templates select="primCategory" />
  </xsl:template>
  <!-- Do nothing for remainder -->
  
  <xsl:template match="version">
    <h2 align="center">
      GLAST Primary Constants, Version 
      <xsl:value-of select="@major" />.<xsl:value-of select="@minor" />
    </h2> 
  </xsl:template>

<xsl:template match="primCategory">
<h2>Primary Category 
  <font color="#f08000">
    <a name="{@name}">
      <xsl:value-of select="@name" />
    </a>
  </font> 
</h2>
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

<xsl:template match="derived" >
    <xsl:apply-templates select="derCategory" />
</xsl:template>

<xsl:template match="derCategory" >
<h2>Derived category
<font color="#f08000">
  <a name="{@name}">
    <xsl:value-of select="@name" />
  </a>
</font> </h2>
   <xsl:call-template name="overview" />
     <table border="1" cellpadding="3">
     <tr align="left" bgcolor="#c0coff">
         <th>Name</th><th>Value</th><th>Description</th>
     </tr>
      <xsl:apply-templates select="const" />
     </table>
<br /><xsl:text> &#xa0;</xsl:text><br />
</xsl:template>

<xsl:template match="const" >
  <tr>
  <td><b><xsl:value-of select="@name" /></b></td>

  <td><xsl:value-of select="@value" />
    <xsl:if test = " @uType = 'length' ">
      <xsl:text> &#xa0;</xsl:text>
      <xsl:value-of select="@unitLength" />
    </xsl:if>
</td>
  <td><xsl:value-of select="notes" /></td>
  </tr>
</xsl:template>

<!--  <xsl:template match="derCategory/overview" ></xsl:template> -->
</xsl:stylesheet>








