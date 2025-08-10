--------------------------------------------------------------------------------
-- This script uses a dedicated mesh generation routine to create artificial  --
-- parameterized spine and spine ER morphologies.                             --
-- It then performs a calcium simulation on the created geometry.             --
-- This script has been used to do the simulations published in:              --
-- Breit et al.: "Spine-to-dendrite calcium modeling discloses relevance for  --
--                precise positioning of ryanodine receptor-containing spine  --
--                endoplasmic reticulum", Scientific Reports (2018)           --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date:   2016-06-22                                                         --
--------------------------------------------------------------------------------

-- for profiler output
--SetOutputProfileStats(true)

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"neuro_collection", "Parmetis"})

-- dimension
dim = 3

-- choose dimension and algebra
InitUG(dim, AlgebraType("CPU", 1))
 
-- choose outfile directory
outPath = util.GetParam("-outName", "output")
outPath = outPath.."/"

-- choice of grid
gridName = util.GetParam("-grid", outPath.."grid/spine.ugx")

-- app length
neckRad = util.GetParamNumber("-neckRad", 0.08)
neckLen = util.GetParamNumber("-neckLen", 0.7)
headRad = util.GetParamNumber("-headRad", 0.29)
headLen = util.GetParamNumber("-headLen", 0.58)
appRad = util.GetParamNumber("-appRad", 0.025)
appLen = util.GetParamNumber("-appLen", 1.175)
appHeadLen = util.GetParamNumber("-appHeadLen", 0.3)
appHeadRad = util.GetParamNumber("-appHeadRad", 0.15)
buildApp = true
if appLen == 0 then
	buildApp = false
end
buildAppHead = false
if appHeadLen > 0 and appHeadRad > 0 then
	buildAppHead = true
end

if ProcRank() == 0 then
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
end

PclDebugBarrierAll()


-- which ER mechanisms are to be activated?
setting = util.GetParam("-setting", "all")
setting = string.lower(setting)
validSettings = {}
validSettings["all"] = 0
validSettings["none"] = 0
validSettings["ip3r"] = 0
validSettings["ryr"] = 0
if (validSettings[setting] == nil) then
    error("Unknown setting " .. setting)
end

withIP3R = true
withRyR = true
withSERCAandLeak = true

if setting == "none" then 
	withIP3R = false
	withRyR = false
	withSERCAandLeak = false
end

if setting == "ip3r" then
	withRyR = false
end

if setting == "ryr" then
	withIP3R = false
end


-- total refinements
numRefs = util.GetParamNumber("-numRefs", 0)

-- choose length of maximal time step during the whole simulation
timeStep = util.GetParamNumber("-tstep", 0.01)

-- choose length of time step at the beginning
-- if not timeStepStart = 2^(-n)*timeStep, take nearest lower number of that form
timeStepStart = util.GetParamNumber("-tstepStart", timeStep)
function log2(x)
	return math.log(x)/math.log(2)
end
startLv =  math.ceil(log2(timeStep/timeStepStart))
timeStepStartNew = timeStep / math.pow(2, startLv)
if (math.abs(timeStepStartNew-timeStepStart)/timeStepStart > 1e-5) then 
	print("timeStepStart argument ("..timeStepStart..") was not admissible; taking "..timeStepStartNew.." instead.")
end
timeStepStart = timeStepStartNew
	
-- choose end time
endTime = util.GetParamNumber("-endTime")
if (endTime == nil)
then
	-- choose number of time steps
	nTimeSteps = util.GetParamNumber("-nTimeSteps", 1)
	endTime = nTimeSteps*timeStep
end

-- choose plotting interval
plotStep = util.GetParamNumber("-pstep", 0.01)

-- choose solver setup
solverID = util.GetParam("-solver", "GS")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0;
validSolverIDs["GS"] = 0;
validSolverIDs["ILU"] = 0;
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")

-- specify -vtk to generate vtk output
verbose = util.HasParamOption("-verbose")

---------------
-- constants --
---------------
-- total cytosolic calbindin concentration
-- (four times the real value in order to simulate four binding sites in one)
totalClb = 4*40.0e-6

-- diffusion coefficients
D_cac = 220.0
D_cae = 220.0
D_ip3 = 280.0
D_clb = 20.0

-- calbindin binding rates
k_bind_clb = 	27.0e06
k_unbind_clb = 	19

-- initial concentrations
ca_ext = 1e-3
ca_cyt_init = 5.0e-08 --4.0e-8
ca_er_init = 2.5e-4
ip3_init = 4.0e-8
clb_init = totalClb / (k_bind_clb/k_unbind_clb*ca_cyt_init + 1)

-- reaction reate IP3
reactionRateIP3 = 0.11

-- equilibrium concentration IP3
equilibriumIP3 = 4.0e-08

-- reation term IP3
reactionTermIP3 = -reactionRateIP3 * equilibriumIP3

IP3Rdensity = 17.3
RYRdensity = 6.0 --3.0 --0.86
leakERconstant = 3.8e-17

-- this is a little bit more complicated, since it must be ensured that
-- the net flux for equilibrium concentrations is zero
-- MUST be adapted whenever any parameterization of ER flux mechanisms is changed!
local v_s = 6.5e-27						-- V_S param of SERCA pump
local k_s = 1.8e-7						-- K_S param of SERCA pump
local j_ip3r = 3.7606194166520605e-23   -- single channel IP3R flux (mol/s) - to be determined via gdb
local j_ryr = 1.1201015633466695e-21    -- single channel RyR flux (mol/s) - to be determined via gdb
				  							-- ryr1: 1.1204582669024472e-21	
local j_leak = ca_er_init-ca_cyt_init	-- leak proportionality factor
	
SERCAdensity = leakERconstant * j_leak
if withIP3R then 
	SERCAdensity = SERCAdensity + IP3Rdensity * j_ip3r
end
if withRyR then
	SERCAdensity = SERCAdensity + RYRdensity * j_ryr
end
SERCAdensity = SERCAdensity / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)

pmcaDensity = 500.0
ncxDensity  = 15.0
--vgccDensity = 1.0

leakPMconstant =  pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s)
				+ ncxDensity *  6.7567567567567566e-23	-- single pump NCX flux (mol/s)
				--+ vgccDensity * (-1.5752042094823713e-25)    -- single channel VGCC flux (mol/s)
				-- *1.5 // * 0.5 for L-type // T-type
if leakPMconstant < 0 then error("PM leak flux is outward for these density settings!") end


-- firing pattern of the synapse
synSubset = 6
caEntryDuration = 0.001 -- was: 0.01


--[[
-- AMPAR
g_max = 0.1
tau = 0.0058

-- NMDAR
g_max = 0.023 -- 0.23 mol/(um^2 s) * (um^3/dm^3)
              -- = 1% of 80pA (synaptic max current) per 0.0181um^2 synapse area (times 1/(2F))
tau = 0.15    -- 150 ms (both values: Destexhe, Mainen & Sejnowski, 1995)
--]]

--arbitrary
g_max = 0.005
tau = 0.01

-- burst of calcium influx for active synapses (~1200 ions)
freq = 50	-- spike train frequency (Hz) (the ineq. 1/freq > caEntryDuration must hold)
nSpikes = 1	-- number of spikes	
function synCurrentDensityCa(x, y, z, t, si)	
	---[[
	-- spike train
	if (si == synSubset and t <= caEntryDuration + (nSpikes - 1) * 1.0/freq) then
        t = t % (1.0/freq)
	end --> now, treat like single spike
	
	-- single spike
	if 	(si == synSubset and t <= caEntryDuration)
	then influx = 0.1 * (1.0 - t/caEntryDuration) ---2e-4 -- accounting for minuscule synapse area
	else influx = 0.0
	end
    return influx
    --]]
    --[[
    if si ~= synSubset then
    	return 0.0
    end
    
    return g_max * math.exp(-t/tau)
    --]]
end

-- burst of ip3 at active synapse (triangular, immediate)
ip3EntryDelay = 0.000
ip3EntryDuration = 0.2
function synCurrentDensityIP3(x, y, z, t, si)
	if (si == synSubset and t > ip3EntryDelay and t <= ip3EntryDelay+ip3EntryDuration)
	then influx = 5e-3 * (1.0 - t/ip3EntryDuration)
	else influx = 0.0
	end
    return influx
end

-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
print("create, refine and distribute domain")

if buildApp then
	reqSubsets = {"cyt", "er", "app", "mem_cyt", "mem_er", "mem_app", "syn", "measZone_dend", "measZone_neck", "measZone_head"}
else
	reqSubsets = {"cyt", "er", "mem_cyt", "mem_er", "syn", "measZone_dend", "measZone_neck", "measZone_head"}
end	
dom = util.CreateDomain(gridName, 0, reqSubsets)

-- in parallel environments: use a load balancer to distribute the grid
-- actual refinement and load balancing after setup of disc.
balancer.partitioner = "parmetis"
balancer.staticProcHierarchy = true
balancer.firstDistLvl = -1
balancer.redistSteps = 0

loadBalancer = balancer.CreateLoadBalancer(dom)
if loadBalancer ~= nil then
	if balancer.partitioner == "parmetis" then
		mu = ManifoldUnificator(dom)
		if buildApp then
			mu:add_protectable_subsets("mem_er, mem_app")
		else
			mu:add_protectable_subsets("mem_er")
		end
		cdgm = ClusteredDualGraphManager()
		cdgm:add_unificator(SiblingUnificator())
		cdgm:add_unificator(mu)
		balancer.defaultPartitioner:set_dual_graph_manager(cdgm)
	end
	balancer.Rebalance(dom, loadBalancer)
end

if numRefs > 0 then	
	refiner = GlobalDomainRefiner(dom)
	for i = 1, numRefs do
		refiner:refine()
	end
end

if loadBalancer ~= nil then
	print("Edge cut on base level: "..balancer.defaultPartitioner:edge_cut_on_lvl(0))
	loadBalancer:estimate_distribution_quality()
	loadBalancer:print_quality_records()
end
print(dom:domain_info():to_string())


--[[
--print("Saving domain grid and hierarchy.")
--SaveDomain(dom, "refined_grid_p" .. ProcRank() .. ".ugx")
--SaveGridHierarchyTransformed(dom:grid(), "refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 2.0)
print("Saving parallel grid layout")
SaveParallelGridLayout(dom:grid(), "parallel_grid_layout_p"..ProcRank()..".ugx", 2.0)
--]]

-- create approximation space
approxSpace = ApproximationSpace(dom)

cytVol = "cyt"
measZones = "measZone_head, measZone_neck, measZone_dend"
cytVol = cytVol .. ", " .. measZones

if buildApp then
	erVol = "er, app"
	erMem = "mem_er, mem_app"
	erMemVec = {"mem_er", "mem_app"}
else
	erVol = "er"
	erMem = "mem_er"
	erMemVec = {"mem_er"}
end
plMem = "mem_cyt, syn"
plMem_vec = {"mem_cyt", "syn"}


outerDomain = cytVol .. ", " .. plMem .. ", " .. erMem
innerDomain = erVol .. ", " .. erMem

approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain)
approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)
approxSpace:add_fct("clb", "Lagrange", 1, outerDomain)

approxSpace:init_top_surface()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()


--------------------------
-- setup discretization --
--------------------------
-- diffusion --
elemDiscCYT = ConvectionDiffusion("ca_cyt", cytVol, "fv1")
elemDiscCYT:set_diffusion(D_cac)

elemDiscER = ConvectionDiffusion("ca_er", erVol, "fv1") 
elemDiscER:set_diffusion(D_cae)

elemDiscIP3 = ConvectionDiffusion("ip3", cytVol, "fv1")
elemDiscIP3:set_diffusion(D_ip3)
elemDiscIP3:set_reaction_rate(reactionRateIP3)
elemDiscIP3:set_reaction(reactionTermIP3)

elemDiscClb = ConvectionDiffusion("clb", cytVol, "fv1")
elemDiscClb:set_diffusion(D_clb)


-- buffering --
elemDiscBuffering = BufferFV1(cytVol)	-- where buffering occurs
elemDiscBuffering:add_reaction(
	"clb",						    -- the buffering substance
	"ca_cyt",						-- the buffered substance
	totalClb,						-- total amount of buffer
	k_bind_clb,					    -- binding rate constant
	k_unbind_clb)				    -- unbinding rate constant


-- er membrane transport systems

-- We pass the function needed to evaluate the flux function here.
-- The order, in which the discrete fcts are passed, is crucial!
ip3r = IP3R({"ca_cyt", "ca_er", "ip3"})
ip3r:set_scale_inputs({1e3,1e3,1e3})
ip3r:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

--ryr = RyR({"ca_cyt", "ca_er"})
ryr = RyRinstat({"ca_cyt", "ca_er"}, erMemVec, approxSpace)
ryr:set_scale_inputs({1e3,1e3})
ryr:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

serca = SERCA({"ca_cyt", "ca_er"})
serca:set_scale_inputs({1e3,1e3})
serca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakER = Leak({"ca_er", "ca_cyt"})
leakER:set_scale_inputs({1e3,1e3})
leakER:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)


discIP3R = MembraneTransportFV1(erMem, ip3r)
discIP3R:set_density_function(IP3Rdensity)

discRyR = MembraneTransportFV1(erMem, ryr)
discRyR:set_density_function(RYRdensity)

discSERCA = MembraneTransportFV1(erMem, serca)
discSERCA:set_density_function(SERCAdensity)

discERLeak = MembraneTransportFV1(erMem, leakER)
discERLeak:set_density_function(1e12*leakERconstant/(1e3)) -- from mol/(um^2 s M) to m/s


-- plasma membrane transport systems
pmca = PMCA({"ca_cyt", ""})
pmca:set_constant(1, ca_ext)
pmca:set_scale_inputs({1e3,1e3})
pmca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ncx = NCX({"ca_cyt", ""})
ncx:set_constant(1, ca_ext)
ncx:set_scale_inputs({1e3,1e3})
ncx:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakPM = Leak({"", "ca_cyt"})
leakPM:set_constant(0, ca_ext)
leakPM:set_scale_inputs({1e3,1e3})
leakPM:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)

--[[
vdcc = VDCC_BG_VM2UG({"ca_cyt", ""}, plMem_vec, approxSpace,
					 "neuronRes/timestep".."_order".. 0 .."_jump"..string.format("%1.1f", 5.0).."_",
					 "%.3f", ".dat", false)
vdcc:set_constant(1, 1.5)
vdcc:set_scale_inputs({1e3,1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_L() --default, but to be sure
vdcc:set_file_times(0.001, 0.0)
vdcc:init(0.0)

vdcc = VDCC_BG_UserData({"ca_cyt", ""}, plMem_vec, approxSpace)
vdcc:set_constant(1, 1.5)
vdcc:set_scale_inputs({1e3,1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_L() --default, but to be sure
vdcc:set_potential_function(-0.065)
vdcc:init(0.0)
--]]

discPMCA = MembraneTransportFV1(plMem, pmca)
discPMCA:set_density_function(pmcaDensity)

discNCX = MembraneTransportFV1(plMem, ncx)
discNCX:set_density_function(ncxDensity)

discPMLeak = MembraneTransportFV1(plMem, leakPM)
discPMLeak:set_density_function(1e9*leakPMconstant / (ca_ext - ca_cyt_init))

--discVDCC = MembraneTransportFV1(plMem, vdcc)
--discVDCC:set_density_function(vdccDensity)


-- synaptic activity
synapseInfluxCa = UserFluxBoundaryFV1("ca_cyt", plMem)
synapseInfluxCa:set_flux_function("synCurrentDensityCa")

synapseInfluxIP3 = UserFluxBoundaryFV1("ip3", plMem)
synapseInfluxIP3:set_flux_function("synCurrentDensityIP3")


------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

domainDisc:add(elemDiscER)
domainDisc:add(elemDiscCYT)
domainDisc:add(elemDiscIP3)
domainDisc:add(elemDiscClb)

domainDisc:add(elemDiscBuffering)

if withIP3R then
	domainDisc:add(discIP3R)
end
if withRyR then
	domainDisc:add(discRyR)
end
if withSERCAandLeak then
	domainDisc:add(discSERCA)
	domainDisc:add(discERLeak)
end

domainDisc:add(discPMCA)
domainDisc:add(discNCX)
domainDisc:add(discPMLeak)
--domainDisc:add(discVDCC)

domainDisc:add(synapseInfluxCa)
domainDisc:add(synapseInfluxIP3)

-- setup time discretization --
timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()

------------------
-- solver setup --
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_base_dir(outPath)
dbgWriter:set_vtk_output(false)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-8)
convCheck:set_verbose(verbose)

if (solverID == "ILU") then
    bcgs_steps = 10000
    ilu = ILU()
    ilu:set_sort(true)
    bcgs_precond = ilu
elseif (solverID == "GS") then
    bcgs_steps = 10000
    bcgs_precond = GaussSeidel()
else -- (solverID == "GMG")
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(timeDisc)
	gmg:set_base_level(0)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	
	-- treat SuperLU problems with Dirichlet constraints by using constrained version
	gmg:set_base_solver(SuperLU())
	
	ilu_gmg = ILU()
	ilu_gmg:set_sort(true)		-- <-- SUPER-important!
	gmg:set_smoother(ilu_gmg)
	gmg:set_smooth_on_surface_rim(true)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	--gmg:set_rap(true) -- causes error in base solver!!
	--gmg:set_debug(GridFunctionDebugWriter(approxSpace))
	
    bcgs_steps = 1000
	bcgs_precond = gmg
end

convCheck:set_maximum_steps(bcgs_steps)

bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(bcgs_precond)
bicgstabSolver:set_convergence_check(convCheck)
--bicgstabSolver:set_debug(dbgWriter)

--- non-linear solver ---
-- convergence check
newtonConvCheck = CompositeConvCheck(approxSpace, 10, 1e-18, 1e-12)
--newtonConvCheck:set_component_check("ca_cyt, ca_er, clb, ip3", 1e-18, 1e-12)
newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)
--newtonConvCheck:set_adaptive(true)

-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)

newtonSolver:init(op)

-------------
-- solving --
-------------

-- get grid function
u = GridFunction(approxSpace)

-- set initial value
InterpolateInner(ca_cyt_init, u, "ca_cyt", 0.0)
InterpolateInner(ca_er_init, u, "ca_er", 0.0)
InterpolateInner(clb_init, u, "clb", 0.0)
InterpolateInner(ip3_init, u, "ip3", 0.0)


-- timestep in seconds
dt = timeStepStart
time = 0.0
step = 0

if (generateVTKoutput) then
	out = VTKOutput()
	out:print(outPath .. "vtk/result", u, step, time)
end

take_measurement(u, time, measZones, "ca_cyt, ip3, clb", outPath .. "meas/data")


-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

if buildApp then
	compute_volume(approxSpace, "measZone_head, measZone_neck, app, syn, mem_er")
else
	compute_volume(approxSpace, "measZone_head, measZone_neck, syn, mem_er")
end

min_dt = timeStep / math.pow(2,15)
cb_interval = 10
lv = startLv
levelUpDelay = 2*caEntryDuration + (nSpikes - 1) * 1.0/freq;
cb_counter = {}
for i=0,startLv do cb_counter[i]=0 end
while endTime-time > 0.001*dt do
	print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
	
	-- setup time Disc for old solutions and timestep
	timeDisc:prepare_step(solTimeSeries, dt)
	
	-- apply newton solver
	if newtonSolver:apply(u) == false
	then
		-- in case of failure:
		print ("Newton solver failed at point in time " .. time .. " with time step " .. dt)
		
		dt = dt/2
		lv = lv + 1
		VecScaleAssign(u, 1.0, solTimeSeries:latest())
		
		-- halve time step and try again unless time step below minimum
		if dt < min_dt
		then 
			print ("Time step below minimum. Aborting. Failed at point in time " .. time .. ".")
			time = endTime
		else
			print ("Trying with half the time step...")
			cb_counter[lv] = 0
		end
	else
		-- update new time
		time = solTimeSeries:time(0) + dt
		
		-- update check-back counter and if applicable, reset dt
		cb_counter[lv] = cb_counter[lv] + 1
		while cb_counter[lv] % (2*cb_interval) == 0 and lv > 0 and (time >= levelUpDelay or lv > startLv) do
			print ("Doubling time due to continuing convergence; now: " .. 2*dt)
			dt = 2*dt;
			lv = lv - 1
			cb_counter[lv] = cb_counter[lv] + cb_counter[lv+1] / 2
			cb_counter[lv+1] = 0
		end
		
		-- plot solution every plotStep seconds
		if (generateVTKoutput) then
			if math.abs(time/plotStep - math.floor(time/plotStep+0.5)) < 1e-5 then
				out:print(outPath .. "vtk/result", u, math.floor(time/plotStep+0.5), time)
				out:write_time_pvd(outPath .. "vtk/result", u)
			end
		end
		
		take_measurement(u, time, measZones, "ca_cyt, ip3, clb", outPath .. "meas/data")
		
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()
		
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is popped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		print("++++++ POINT IN TIME  " .. math.floor(time/dt+0.5)*dt .. "s  END ++++++++");
	end

end

