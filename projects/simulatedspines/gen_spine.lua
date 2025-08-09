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
out_name = util.GetParam("-outName","output")
file_exists = hlp.dir_exists(out_name);
if not file_exists then
    print("Output directory exists: "..tostring(file_exists))
    os.execute("mkdir -p " .. out_name)
    print("Made directory: "..out_name)
else
    print("Directory: "..out_name.." exists!")
end


