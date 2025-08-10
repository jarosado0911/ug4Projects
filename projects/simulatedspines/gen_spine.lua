--------------------------------------------------------------------------------
-- This script uses a dedicated mesh generation routine to create artificial  --
-- parameterized spine and spine ER morphologies.                             --
-- I have borrowed and modify from an original script Markus Breit designed   --
-- The details of the script where published in the article below             --
-- Breit et al.: "Spine-to-dendrite calcium modeling discloses relevance for  --
--                precise positioning of ryanodine receptor-containing spine  --
--                endoplasmic reticulum", Scientific Reports (2018)           --
--                                                                            --
-- Author: James Rosado                                                       --
-- Date:   2025-09-25                                                         --
--------------------------------------------------------------------------------

-- this is a helper module I made
local hlp = require "helperf"

-- this will initialize the lua script
hlp.initialize_script()

-- choose outfile directory
outName = util.GetParam("-outName","output")
file_exists = hlp.dir_exists(outName);
if not file_exists then
    print("Output directory exists: "..tostring(file_exists))
    os.execute("mkdir -p " .. outName)
    print("Made directory: "..outName)
else
    print("Directory: "..outName.." exists!")
end

-- set the output grid file name
gridName = outName .. "/artificalspine.ugx"

-- now let us collect some parameters from the user
-- app length
neckRad = util.GetParamNumber("-neckRad", 0.08)
neckLen = util.GetParamNumber("-neckLen", 0.7)
headRad = util.GetParamNumber("-headRad", 0.29)
headLen = util.GetParamNumber("-headLen", 0.58)
appRad = util.GetParamNumber("-appRad", 0.025)
appLen = util.GetParamNumber("-appLen", 0.5)
appHeadLen = util.GetParamNumber("-appHeadLen", 0.0)
appHeadRad = util.GetParamNumber("-appHeadRad", 0.0)
buildApp = true
if appLen == 0 then
	buildApp = false
end
buildAppHead = false
if appHeadLen > 0 and appHeadRad > 0 then
	buildAppHead = true
end

-- print a summary of parameters
print("\n  Summary of Morphology Specifications")
print("------------------------------------------")
print("\tneckRad:      "..tostring(neckRad))
print("\tneckLen:      "..tostring(neckLen))
print("\theadRad:      "..tostring(headRad))
print("\theadLen:      "..tostring(headLen))
print("\tappRad:       "..tostring(appRad))
print("\tappLen:       "..tostring(appLen))
print("\tappHeadLen:   "..tostring(appHeadLen))
print("\tappHeadRad:   "..tostring(appHeadRad))
print("\tbuildApp:     "..tostring(buildApp))
print("\tbuildAppHead: "..tostring(buildAppHead))

-- you will need tetgen to run this method
BuildDendrite(
	{
		0.45,	-- cytosol radius (um)
		0.11,	-- ER radius (um)
		10.0,	-- dendrite length (um)
		5.0,	-- spine position (um)
		appRad,	-- app neck radius (um)
		appLen,	-- app neck length (um)
		appHeadRad - appRad,	-- app head radius (in addition to app neck radius) (um)
		appHeadLen,	-- app head length (um)
		neckRad,	-- spine neck radius (um)
		neckLen,-- spine neck length (um)
		headRad - neckRad,-- spine head radius (in addition to spine neck radius) (um)
		headLen	-- spine head length (um)
	},
	{
		true,	-- synapse?
		true,	-- ER?
		buildApp,	-- app?
		false,	-- never use that
		buildAppHead	-- app head?
	},
	gridName
)
print("Building the dendrite")
