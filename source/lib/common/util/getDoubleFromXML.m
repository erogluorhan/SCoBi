function output = getDoubleFromXML( xmlFile, varName )

output = str2double( xmlFile.getElementsByTagName(varName).item(0).getFirstChild.getData );

end