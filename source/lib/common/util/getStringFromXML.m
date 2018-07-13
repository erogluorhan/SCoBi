function output = getStringFromXML( xmlFile, varName )

output = char( xmlFile.getElementsByTagName(varName).item(0).getFirstChild.getData );

end