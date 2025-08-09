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

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")

AssertPluginsLoaded({"neuro_collection"})