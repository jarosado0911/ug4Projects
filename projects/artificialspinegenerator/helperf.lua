local helperf = {}

function helperf.file_exists(file)
  -- some error codes:
  -- 13 : EACCES - Permission denied
  -- 17 : EEXIST - File exists
  -- 20	: ENOTDIR - Not a directory
  -- 21	: EISDIR - Is a directory
  --
  local isok, errstr, errcode = os.rename(file, file)
  if isok == nil then
     if errcode == 13 then 
        -- Permission denied, but it exists
        return true
     end
     return false
  end
  return true
end
 
function helperf.dir_exists(path)
  return helperf.file_exists(path .. "/")
end

function helperf.initialize_script()
	-- initialize the script
	print("Loading utilities...")
	ug_load_script("ug_util.lua")
	print("Loading ug_util.lua ... success!")
	print("Asserting ProMesh plugin...")
	AssertPluginsLoaded({"neuro_collection","ProMesh"})
	--print("Asserting ProMesh plugin ... success!")

	-- choose dimension and algebra
	InitUG(3, AlgebraType("CPU", 1))
end

-- must have this last line below!
return helperf