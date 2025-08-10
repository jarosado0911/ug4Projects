--------------------------------------------------------------------------------
-- Example script for calcium simulations on a 3d reconstructed spine without --
-- spine endoplasmic reticulum using a simple adaptive time stepping method.  --
--                                                                            --
-- Author: Markus Breit                                                       --
-- Date: 2019-05-17                                                           --
--------------------------------------------------------------------------------

-- for profiler output
SetOutputProfileStats(false)

-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"neuro_collection", "MembranePotentialMapping"})

-- choose dimension and algebra
InitUG(3, AlgebraType("CPU", 1))

-- choice of grid
gridName = util.GetParam("-grid", "calciumDynamics_app/grids/reconstructed_spine.ugx")

-- total refinements
numRefs = util.GetParamNumber("-numRefs", 0)

-- choose length of maximal time step during the whole simulation
timeStep = util.GetParamNumber("-tstep", 1e-04)

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
if endTime == nil
then
	-- choose number of time steps
	nTimeSteps = util.GetParamNumber("-nTimeSteps", 1)
	endTime = nTimeSteps*timeStep
end

-- choose solver setup
solverID = util.GetParam("-solver", "GS")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0
validSolverIDs["GS"] = 0
validSolverIDs["ILU"] = 0
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- specify -verbose to display linear solver output
verbose = util.HasParamOption("-verbose")
 
-- choose outfile directory
fileName = util.GetParam("-outName", "CD/reconstructed")
fileName = fileName.."/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")

-- choose plotting interval
plotStep = util.GetParamNumber("-pstep", timeStep)

---------------
-- constants --
---------------
-- total cytosolic calbindin concentration
-- (four times the real value in order to simulate four binding sites in one)
totalClb = 4*40.0e-6

-- diffusion coefficients
D_cac = 220.0
D_clb = 20.0

-- calbindin binding rates
k_bind_clb = 	27.0e06
k_unbind_clb = 	19

-- initial concentrations
ca_cyt_init = 5.0e-08
ca_ext = 1.0e-3
clb_init = totalClb / (k_bind_clb/k_unbind_clb*ca_cyt_init + 1)

pmcaDensity = 500.0
ncxDensity  = 15.0
vdccDensity = 1.0

leakPMconstant =  pmcaDensity * 6.9672131147540994e-24	-- single pump PMCA flux (mol/s) at 1mM ext. Ca2+
				+ ncxDensity *  6.7567567567567566e-23	-- single pump NCX flux (mol/s) at 1mM ext. Ca2+
				- vdccDensity * 1.435114751437757757e-28  -- single channel VGCC flux (mol/s) at 1mM ext. Ca2+ and V_m = -0.07V
				-- *1.5 // * 0.5 for L-type // T-type
if leakPMconstant < 0 then error("PM leak flux is outward for these density settings!") end


-- firing pattern of the synapse
synStartTime = 0.0
caEntryDuration = 0.01

-- burst of calcium influx for active synapses (~1200 ions)
freq = 50      -- spike train frequency (Hz) (the ineq. 1/freq > caEntryDuration must hold)
nSpikes = 1    -- number of spikes	
function neumannBndCa(x, y, z, t, si)
	if t < synStartTime then
		return 0.0
	end
	
	-- spike train
	if t <= synStartTime + caEntryDuration + (nSpikes - 1) * 1.0/freq then
        t = t % (1.0/freq)
	end --> now, treat like single spike
	
	-- single spike
	if 	t <= caEntryDuration then
		influx = 0.0086*(1.0 - t/caEntryDuration)
	else
		influx = 0.0
	end
	
    return influx
end


-------------------------------
-- setup approximation space --
-------------------------------
-- create, load, refine and distribute domain
print("create, refine and distribute domain")

reqSubsets = {"cyt", "pm", "syn"} --, "meas_dend", "meas_neck", "meas_head"}
dom = util.CreateDomain(gridName, 0, reqSubsets)

print(dom:domain_info():to_string())

-- create approximation space
approxSpace = ApproximationSpace(dom)

cytVol = "cyt"
--measZones = "meas_head, meas_neck, meas_dend"
measZones = "cyt"
--cytVol = cytVol .. ", " .. measZones

plMem = "pm, syn"
plMem_vec = {"pm", "syn"}

outerDomain = cytVol .. ", " .. plMem

approxSpace:add_fct("ca_cyt", "Lagrange", 1)
approxSpace:add_fct("clb", "Lagrange", 1)
approxSpace:add_fct("m", "Lagrange", 1, plMem)

approxSpace:init_levels()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

--------------------------
-- setup discretization --
--------------------------
-- diffusion --
elemDiscCYT = ConvectionDiffusion("ca_cyt", cytVol, "fv1")
elemDiscCYT:set_diffusion(D_cac)

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
leakPM:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)


exper_type = 0 -- for VDCC only
--exper_type = 1 -- VDCC and Ca-influx
--exper_type = 2 -- VDCC and AMPAR (vdcc on pm only)
--exper_type = 3 -- VDCC and NMDAR (vdcc on pm only)

vdccMode = 0  -- 0 for script function voltage, 1 for data file voltage
if vdccMode == 0 then
	apSignal = ActionPotentialTrain(0.0, 0.02, 50, -70.0)
	function membranePotential(x, y, z, t, si)
		-- 80% of the intensity of the real AP
		-- return 1e-3*(-70.0 + 0.8*(apSignal:membrane_potential(t) + 70.0))
		if exper_type == 0 then
			return 1e-3*(70)
		elseif exper_type == 1 then
			return 1e-3*(-70)
		elseif exper_type == 2 then
			return 1e-3*(-70)
		elseif exper_type == 3 then
			return 1e-3*(-70)
		end

	end
	--limit vdcc to pm when AMPAR on
	if exper_type == 2 then
		vdcc = VDCC_BG_UserData({"ca_cyt", "", "m"}, {"pm"}, approxSpace)
		vdcc:set_potential_function("membranePotential")
	--limit vdcc to pm when NMDAR on	
	elseif exper_type == 3 then
		vdcc = VDCC_BG_UserData({"ca_cyt", "", "m"}, {"pm"}, approxSpace)
		vdcc:set_potential_function("membranePotential")
	--else put vdcc everywhere
	else
		vdcc = VDCC_BG_UserData({"ca_cyt", "", "m"}, plMem_vec, approxSpace)
		vdcc:set_potential_function("membranePotential")
	end
elseif vdccMode == 1 then
	--limit vdcc to pm when AMPAR on
	if exper_type == 2 then
		vdcc = VDCC_BG_VM2UG({"ca_cyt", "", "m"}, {"pm"}, approxSpace,
				"voltageData/vm_", "%0.5f", ".dat", false)
		vdcc:set_file_times(1e-5, 0.0) -- file interval is 0.01ms, starting at 0.0
	--limit vdcc to pm when NMDAR on
	elseif exper_type == 3 then
		vdcc = VDCC_BG_VM2UG({"ca_cyt", "", "m"}, {"pm"}, approxSpace,
				"voltageData/vm_", "%0.5f", ".dat", false)
		vdcc:set_file_times(1e-5, 0.0) -- file interval is 0.01ms, starting at 0.0
	--else put vdcc everywhere
	else
		vdcc = VDCC_BG_VM2UG({"ca_cyt", "", "m"}, plMem_vec, approxSpace,
				"voltageData/vm_", "%0.5f", ".dat", false)
		vdcc:set_file_times(1e-5, 0.0) -- file interval is 0.01ms, starting at 0.0
	end
end
vdcc:set_constant(1, ca_ext)
vdcc:set_scale_inputs({1e3, 1e3, 1.0})
vdcc:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
vdcc:set_channel_type_L()
vdcc:init(0.0)

ampaMode = 0  -- 0 for script function voltage, 1 for data file voltage
if ampaMode == 0 then
	ampa_apSignal = ActionPotentialTrain(0.0, 0.02, 50, -70.0)
	function ampa_membranePotential(x, y, z, t, si)
		-- 80% of the intensity of the real AP
		-- return 1e-3*(-70.0 + 0.8*(apSignal:membrane_potential(t) + 70.0))
		return 1e-3*(70)
	end
	
	ampa = VDCC_BG_UserData({"ca_cyt", "", "m"}, {"syn"}, approxSpace)
	ampa:set_potential_function("ampa_membranePotential")
elseif ampaMode == 1 then
	ampa = VDCC_BG_VM2UG({"ca_cyt", "", "m"}, {"syn"}, approxSpace,
			"voltageData/vm_", "%0.5f", ".dat", false)
	ampa:set_file_times(1e-5, 0.0) -- file interval is 0.01ms, starting at 0.0
end
ampa:set_constant(1, ca_ext)
ampa:set_scale_inputs({1e3, 1e3, 1.0})
ampa:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)
ampa:set_channel_type_L()
ampa:init(0.0)

nmdar = NMDAR({"","ca_cyt"})
nmdar:set_constant(0,ca_ext)
nmdar:set_scale_inputs({1e3,1e3})
nmdar:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

discPMCA = MembraneTransportFV1(plMem, pmca)
discPMCA:set_density_function(pmcaDensity)

discNCX = MembraneTransportFV1(plMem, ncx)
discNCX:set_density_function(ncxDensity)

discPMLeak = MembraneTransportFV1(plMem, leakPM)
discPMLeak:set_density_function(leakPMconstant / (1e3*(ca_ext - ca_cyt_init)))

amparDensity = 100.0
discAMPA = MembraneTransportFV1("syn", ampa)
discAMPA:set_density_function(amparDensity)

nmdarDensity = 55.0
discNMDAR = MembraneTransportFV1("syn",nmdar)
discNMDAR:set_density_function(nmdarDensity)

if exper_type == 2 then
	discVDCC = MembraneTransportFV1("pm", vdcc)
	discVDCC:set_density_function(vdccDensity)
elseif exper_type == 3 then
	discVDCC = MembraneTransportFV1("pm",vdcc)
	discVDCC:set_density_function(vdccDensity)
else
	discVDCC = MembraneTransportFV1(plMem, vdcc)
	discVDCC:set_density_function(vdccDensity)
end
-- synaptic activity
synapseInfluxCa = UserFluxBoundaryFV1("ca_cyt", "syn")
synapseInfluxCa:set_flux_function("neumannBndCa")


------------------------------------------
-- setup complete domain discretization --
------------------------------------------
domainDisc = DomainDiscretization(approxSpace)

domainDisc:add(elemDiscCYT)
domainDisc:add(elemDiscClb)

domainDisc:add(elemDiscBuffering)

domainDisc:add(discPMCA)
domainDisc:add(discNCX)
domainDisc:add(discPMLeak)

if exper_type == 0 then
	domainDisc:add(discVDCC)
	domainDisc:add(vdcc)  -- for gating param discretization
elseif exper_type == 1 then
	domainDisc:add(discVDCC)
	domainDisc:add(vdcc)  -- for gating param discretization
	domainDisc:add(synapseInfluxCa)
elseif exper_type == 2 then
	domainDisc:add(discVDCC)
	domainDisc:add(vdcc)  -- for gating param discretization
	domainDisc:add(discAMPA)
	domainDisc:add(ampa)
elseif exper_type == 3 then
	domainDisc:add(discVDCC)
	domainDisc:add(vdcc)  -- for gating param discretization
	domainDisc:add(discNMDAR)
end

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
dbgWriter:set_base_dir(fileName)
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
    --ilu:set_inversion_eps(1.0e-10)
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
	gmg:set_smoother(ilu_gmg)
	gmg:set_smooth_on_surface_rim(true)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	gmg:set_rap(true)
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
newtonConvCheck = CompositeConvCheck(approxSpace, 10, 1e-12, 1e-08)
--newtonConvCheck:set_component_check("ca_cyt, clb", 1e-18, 1e-12)
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
Interpolate(ca_cyt_init, u, "ca_cyt", 0.0)
Interpolate(clb_init, u, "clb", 0.0)
if exper_type == 0 then
	vdcc:calculate_steady_state(u, -0.07)
elseif exper_type == 1 then
	vdcc:calculate_steady_state(u, -0.07)
elseif exper_type == 2 then
	vdcc:calculate_steady_state(u, -0.07)
	ampa:calculate_steady_state(u, -0.07)
elseif exper_type == 3 then
	vdcc:calculate_steady_state(u, -0.07)
end

-- timestep in seconds
dt = timeStepStart
time = 0.0
step = 0

if generateVTKoutput then
	out = VTKOutput()
	out:print(fileName .. "vtk/result", u, step, time)
	if exper_type == 0 then
		--vdcc:export_membrane_potential_to_vtk(fileName .. "vtk_Voltage/vdcc_voltage", step, time)
	elseif exper_type == 2 then
		ampa:export_membrane_potential_to_vtk(fileName .. "vtk_Voltage/ampa_voltage", step, time)
	end
end

take_measurement(u, time, measZones, "ca_cyt, clb", fileName .. "meas/data")

--compute_volume(approxSpace, "meas_head, meas_neck, syn")

-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

min_dt = timeStep / math.pow(2,15)
cb_interval = 10
lv = startLv
levelUpDelay = caEntryDuration
cb_counter = {}
for i=0,startLv do cb_counter[i]=0 end
while endTime-time > 0.001*dt do
	print("++++++ POINT IN TIME  " .. math.floor((time+dt)/dt+0.5)*dt .. "s  BEGIN ++++++")
		
	-- setup time disc for old solutions and time step
	timeDisc:prepare_step(solTimeSeries, dt)
	print("here")
	if newtonSolver:apply(u) == false then
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
		if generateVTKoutput then
			if math.abs(time/plotStep - math.floor(time/plotStep+0.5)) < 1e-5 then
				out:print(fileName .. "vtk/result", u, math.floor(time/plotStep+0.5), time)

				if exper_type == 0 then
					--vdcc:export_membrane_potential_to_vtk(fileName .. "vtk_Voltage/vdcc_voltage", math.floor(time/plotStep+0.5), time)
				elseif exper_type == 2 then
					ampa:export_membrane_potential_to_vtk(fileName .. "vtk_Voltage/ampa_voltage", math.floor(time/plotStep+0.5), time)
				end
			end
		end
		
		take_measurement(u, time, measZones, "ca_cyt, clb", fileName .. "meas/data")
		
		-- get oldest solution
		oldestSol = solTimeSeries:oldest()
		
		-- copy values into oldest solution (we reuse the memory here)
		VecScaleAssign(oldestSol, 1.0, u)
		
		-- push oldest solutions with new values to front, oldest sol pointer is popped from end
		solTimeSeries:push_discard_oldest(oldestSol, time)
		
		print("++++++ POINT IN TIME  " .. math.floor(time/dt+0.5)*dt .. "s  END ++++++++");
	end
end

-- end timeseries, produce gathering file
if generateVTKoutput then
	out:write_time_pvd(fileName .. "vtk/result", u)
end

