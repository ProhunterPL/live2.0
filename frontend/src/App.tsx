import React, { useState, useEffect, useRef } from 'react'
import { Play, Pause, Square, RotateCcw, Database } from 'lucide-react'
import HeatmapCanvas from './components/HeatmapCanvas'
import Controls from './components/Controls'
import GraphPreview from './components/GraphPreview'
import PerformancePanel from './components/PerformancePanel'
import NoveltyPanel from './components/NoveltyPanel'
import APIv1Jobs from './components/APIv1Jobs'
import { SimulationAPI } from './lib/api.ts'
import { WebSocketClient } from './lib/ws.ts'
import type { SimulationData, SimulationStatus, Metrics } from './lib/types.ts'
import logoImage from './assets/logo.jpg'
import headerImage from './assets/header.jpg'

const App: React.FC = () => {
  const [simulationId, setSimulationId] = useState<string | null>(null)
  const [status, setStatus] = useState<SimulationStatus | null>(null)
  const [simulationData, setSimulationData] = useState<SimulationData | null>(null)
  const [metricsSnapshot, setMetricsSnapshot] = useState<Metrics | null>(null)
  const [performanceData, setPerformanceData] = useState<any>(null)
  const [isConnected, setIsConnected] = useState(false)
  const [selectedSubstance, setSelectedSubstance] = useState<string | null>(null)
  const [mode, setMode] = useState<'preset_prebiotic' | 'open_chemistry'>('open_chemistry')
  const [view, setView] = useState<'simulation' | 'api-v1'>('simulation')
  const [availableSpecies, setAvailableSpecies] = useState<string[]>([])
  const [runtimeStartMs, setRuntimeStartMs] = useState<number | null>(null)
  const [runtimeAccumulatedMs, setRuntimeAccumulatedMs] = useState<number>(0)
  const [runtimeNowMs, setRuntimeNowMs] = useState<number>(Date.now())
  
  const api = useRef(new SimulationAPI()).current
  const wsClient = useRef(new WebSocketClient()).current
  
  // FIX: Update runtime clock every second
  useEffect(() => {
    const interval = setInterval(() => {
      setRuntimeNowMs(Date.now())
    }, 1000) // Update every second
    
    return () => clearInterval(interval)
  }, [])
  const wsHandlersRef = useRef<Record<string, Function>>({})

  useEffect(() => {
    if (simulationData?.step_count !== undefined) {
      console.log('✅ Frontend step update:', simulationData.step_count)
      const legendEl = document.querySelector('.legend-column .legend-info strong')
      if (legendEl) {
        legendEl.closest('.legend-column')?.scrollIntoView({ behavior: 'smooth', block: 'nearest' })
      }
    }
  }, [simulationData?.step_count])

  const attachWsListener = (event: string, handler: Function) => {
    const handlers = wsHandlersRef.current
    const existing = handlers[event]
    if (existing) {
      wsClient.off(event, existing)
    }
    handlers[event] = handler
    wsClient.on(event, handler)
  }

  const detachAllWsListeners = () => {
    const handlers = wsHandlersRef.current
    Object.entries(handlers).forEach(([event, handler]) => {
      wsClient.off(event, handler)
    })
    wsHandlersRef.current = {}
  }

  const initializedRef = useRef(false)
  useEffect(() => {
    // NAPRAWIONE: Zawsze łącz się z najnowszą symulacją po odświeżeniu
    // Nie blokujemy przy F5, tylko przy React StrictMode double-mount
    const wasInitialized = initializedRef.current
    initializedRef.current = true
    
    // Po odświeżeniu strony (F5) - zawsze sprawdź current simulation
    if (!wasInitialized) {
      initializeSimulation()
    } else {
      // Jeśli to React StrictMode double-mount, ignoruj
      console.log('Skipping duplicate React StrictMode mount')
    }
    
    return () => {
      // Cleanup: remove listeners and disconnect without stopping simulation
      detachAllWsListeners()
      wsClient.disconnect()
      // Reset na wypadek unmount
      initializedRef.current = false
    }
  }, [wsClient])

  // NOWE: Automatycznie wykryj nową symulację i przełącz się
  useEffect(() => {
    const checkForNewSimulation = async () => {
      try {
        console.log('🔍 Checking for existing simulation...')
        const currentSim = await api.getCurrentSimulation()
        console.log('📊 Current simulation result:', currentSim)
        
        // OPTIMIZATION: Clear old simulation ID if it doesn't exist on backend
        if (simulationId && currentSim.exists && simulationId !== currentSim.simulation_id) {
          console.log(`🧹 Clearing old simulation ID: ${simulationId}`)
          setSimulationId(null)
          setStatus(null)
          setIsConnected(false)
          wsClient.disconnect()
        }
        
        if (currentSim.exists && currentSim.simulation_id) {
          // Jeśli backend ma symulację a my nie, lub ma INNĄ symulację
          if (!simulationId || simulationId !== currentSim.simulation_id) {
            console.log(`🔄 Auto-switching to simulation: ${currentSim.simulation_id}`)
            setSimulationId(currentSim.simulation_id)
            // OPTIMIZATION: Only connect if not already connected to this simulation
            if (!isConnected || simulationId !== currentSim.simulation_id) {
              await connectWebSocket(currentSim.simulation_id)
            }
            setRuntimeAccumulatedMs(0)
            // Sprawdź czy symulacja jest uruchomiona
            const status = await api.getSimulationStatus(currentSim.simulation_id)
            console.log('📊 Simulation status:', status)
            if (!status.is_running) {
              console.log('▶️ Starting existing simulation')
              await startSimulation(currentSim.simulation_id)
              setRuntimeStartMs(Date.now())  // FIX: Start runtime clock
            } else {
              console.log('✅ Simulation already running, starting runtime clock')
              setRuntimeStartMs(Date.now())  // FIX: Start runtime clock for already running simulation
            }
          }
        } else {
          console.log('❌ No existing simulation found')
        }
      } catch (error) {
        console.error('Failed to check for new simulation:', error)
      }
    }
    
    // Sprawdź przy starcie
    checkForNewSimulation()
    
    // OPTIMIZATION: Check every 30s instead of 10s to reduce reconnection frequency
    const interval = setInterval(checkForNewSimulation, 30000)
    return () => clearInterval(interval)
  }, [simulationId, api, isConnected])

  // Poll status periodically so the time advances even if WS data is delayed
  useEffect(() => {
    if (!simulationId) return
    const interval = setInterval(() => {
      updateStatus(simulationId)
    }, 20000)  // Poll status every 20s; WS provides live data
    return () => clearInterval(interval)
  }, [simulationId, api])

  // Runtime clock tick
  useEffect(() => {
    const tick = setInterval(() => setRuntimeNowMs(Date.now()), 200)
    return () => clearInterval(tick)
  }, [])

  // Recreate simulation when mode changes
  useEffect(() => {
    if (!initializedRef.current) return
    const recreate = async () => {
      try {
        // Disconnect current WS
        wsClient.disconnect()
        setIsConnected(false)
        // Reset runtime
        setRuntimeAccumulatedMs(0)
        setRuntimeStartMs(null)
        // Create new simulation in selected mode
        const response = await api.createSimulation({
          config: {  // Note: using 'as any' below to allow extra thermodynamic validation fields
            grid_height: 128,
            grid_width: 128,
            mode,
            max_particles: 5000,
            max_time: 1000,
            dt: 0.005,  // NAPRAWIONE: było 0.035 (7× za dużo!) - stabilność numeryczna
            energy_decay: 0.96,  // Wolniejsze wygaszanie z 0.95
            energy_threshold: 0.1,
            // Energy system parameters - OPTIMIZED FOR MORE DYNAMICS
            pulse_every: 24,  // Częstsze pulsy (co 24 kroki)
            pulse_radius: 32.0,  // Większy promień (32 jednostki)
            pulse_amplitude: 8.0,  // Większa amplituda (8.0)
            diffuse_D: 0.5,  // Szybsza dyfuzja energii
            target_energy: 0.5,  // Wyższe tło energii
            thermostat_alpha: 0.01,  // Silniejszy termostat
            // Thermodynamic validation - SMART VALIDATION (GROMACS/NAMD best practices)
            // Energy+Momentum: every call (~2ms), M-B: every 20k steps, Entropy: every 50k steps
            enable_thermodynamic_validation: true,  // Włączone - smart validation (5,650× szybciej!)
            validate_every_n_steps: 1000,  // Częściej ale z szybkimi testami (Energy+Momentum)
            // Performance optimization parameters
            energy_update_interval: 5,  // Update energy every 5 steps
            metrics_update_interval: 10,  // Update metrics every 10 steps
            diagnostics_frequency: 20,  // Log diagnostics every 20 steps
            // Particle and binding parameters
            particle_radius: 0.5,
            binding_threshold: 0.25,  // NAPRAWIONE: Obniżone z 0.68 dla łatwiejszego tworzenia wiązań
            unbinding_threshold: 0.15,  // NAPRAWIONE: Obniżone z 0.3 dla stabilniejszych wiązań
            // Mutation parameters - OPTIMIZED FOR MORE NOVELTY
            p_mut_base: 5e-4,  // Zwiększone mutacje (5x więcej)
            p_mut_gain: 20.0,  // Zwiększone z 14.0
            attr_sigma: 0.12,  // Silniejsza perturbacja
            // Open chemistry parameters
            vmax: 12.0,  // Zwiększona prędkość maksymalna
            neighbor_radius: 4.0,  // Większy promień sąsiedztwa
            rebuild_neighbors_every: 6,  // Częstsza przebudowa sąsiadów
            clamp_density_per_cell: 64,  // NOWE: maksymalna gęstość na komórkę
            // Other parameters
            novelty_window: 100,
            min_cluster_size: 2,
            vis_frequency: 5,  // Rzadsze renderowanie
            log_frequency: 100,
            seed: Math.floor(Math.random() * 1000000)
          } as any,  // Allow extra config fields for backend (Pydantic accepts additional fields)
          mode
        })
        if (response.success && response.simulation_id) {
          setSimulationId(response.simulation_id)
          await connectWebSocket(response.simulation_id)
          // Auto-start simulation after creation
          setRuntimeAccumulatedMs(0)
          await startSimulation(response.simulation_id)
          // startSimulation already sets runtimeStartMs to Date.now()
        }
      } catch (err) {
        console.error('Failed to recreate simulation:', err)
      }
    }
    recreate()
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [mode])

  const initializeSimulation = async () => {
    try {
      // Check if simulation already exists
      const currentSim = await api.getCurrentSimulation()
      if (currentSim.exists && currentSim.simulation_id) {
        console.log('Using existing simulation:', currentSim.simulation_id)
        setSimulationId(currentSim.simulation_id)
        await connectWebSocket(currentSim.simulation_id)
        setRuntimeAccumulatedMs(0)
        setRuntimeStartMs(null)
        return
      }

      // Create new simulation only if none exists
      const response = await api.createSimulation({
        config: {
          grid_height: 128,
          grid_width: 128,
          mode,
          max_particles: 5000,
          max_time: 1000,
          dt: 0.005,  // NAPRAWIONE: było 0.035 (7× za dużo!) - stabilność numeryczna
          energy_decay: 0.96,  // Wolniejsze wygaszanie z 0.99
          energy_threshold: 0.01,
          pulse_every: 24,  // Częstsze pulsy (co 24 kroki)
          pulse_radius: 32.0,  // Większy promień (32 jednostki)
          pulse_amplitude: 8.0,  // Większa amplituda (8.0)
          diffuse_D: 0.5,  // Szybsza dyfuzja energii
          target_energy: 0.5,  // Wyższe tło energii
          thermostat_alpha: 0.01,  // Silniejszy termostat
          // Performance optimization parameters
          energy_update_interval: 5,  // Update energy every 5 steps
          metrics_update_interval: 10,  // Update metrics every 10 steps
          diagnostics_frequency: 20,  // Log diagnostics every 20 steps
          particle_radius: 0.5,
          binding_threshold: 0.25,  // NAPRAWIONE: Obniżone z 0.68 dla łatwiejszego tworzenia wiązań
          unbinding_threshold: 0.15,  // NAPRAWIONE: Obniżone z 0.3 dla stabilniejszych wiązań
          novelty_window: 100,
          min_cluster_size: 2,
          vis_frequency: 5,  // Rzadsze renderowanie
          log_frequency: 100,
          p_mut_base: 5e-4,  // Zwiększone mutacje (5x więcej)
          p_mut_gain: 20.0,  // Zwiększone z 14.0
          attr_sigma: 0.12,  // Silniejsza perturbacja
          // Open chemistry parameters
          vmax: 12.0,  // Zwiększona prędkość maksymalna
          neighbor_radius: 3.2,  // NOWE: promień sąsiedztwa
          rebuild_neighbors_every: 8,  // NOWE: częstotliwość odbudowy sąsiedztwa
          clamp_density_per_cell: 64,  // NOWE: maksymalna gęstość na komórkę
          seed: Math.floor(Math.random() * 1000000)
        },
        mode
      })

      if (response.success && response.simulation_id) {
        setSimulationId(response.simulation_id)
        await connectWebSocket(response.simulation_id)
        // Auto-start simulation after creation
        setRuntimeAccumulatedMs(0)
        await startSimulation(response.simulation_id)
        // startSimulation already sets runtimeStartMs to Date.now()
      }
    } catch (error) {
      console.error('Failed to initialize simulation:', error)
    }
  }

  const connectWebSocket = async (id: string | null) => {
    try {
      if (!id) {
        console.error('connectWebSocket called without simulationId')
        return
      }
      
      // Check if simulation exists before connecting
      const exists = await api.checkSimulationExists(id)
      if (!exists) {
        console.warn(`Simulation ${id} not found on backend. Clearing simulation ID.`)
        // Don't create new simulation automatically - let user decide
        setSimulationId(null)
        setStatus(null)
        setIsConnected(false)
        return
      }
      
      // Pass only the simulationId; ws client will build the full URL
      detachAllWsListeners()
      wsClient.disconnect()
      await wsClient.connect(id)
      setIsConnected(true)
      
      const handleSimulationData = (data: SimulationData) => {
        setSimulationData(prev => {
          if (!prev || prev.step_count !== data.step_count) {
            console.log('🎯 Updating simulationData to step', data.step_count)
          }
          return data
        })
        if (data.metrics) {
          setMetricsSnapshot(data.metrics)
        }
        
        // Update performance data
        if (data.performance) {
          setPerformanceData(data.performance)
        }
        
        // Extract available species for preset mode
        if (data.concentrations && Object.keys(data.concentrations).length > 0) {
          setAvailableSpecies(Object.keys(data.concentrations))
        }
      }

      const handleDisconnected = () => {
        setIsConnected(false)
      }

      const handleError = (error: Error) => {
        console.error('WebSocket error:', error)
        setIsConnected(false)
      }

      const handleSimulationNotFound = () => {
        console.warn('Simulation not found on backend. Clearing simulation ID.')
        setIsConnected(false)
        setSimulationId(null)
        setStatus(null)
        // Don't create new simulation automatically
      }

      const handleReconnectFailed = () => {
        console.error('Failed to reconnect after multiple attempts')
        setIsConnected(false)
      }

      attachWsListener('simulation_data', handleSimulationData)
      attachWsListener('disconnected', handleDisconnected)
      attachWsListener('error', handleError)
      attachWsListener('simulation_not_found', handleSimulationNotFound)
      attachWsListener('reconnect_failed', handleReconnectFailed)
    } catch (error) {
      console.error('Failed to connect WebSocket:', error)
    }
  }

  const startSimulation = async (id: string) => {
    try {
      await api.startSimulation(id)
      await updateStatus(id)
      await refreshMetrics(id)
      // start wall-clock runtime
      setRuntimeStartMs(Date.now())
    } catch (error) {
      console.error('Failed to start simulation:', error)
    }
  }

  const pauseSimulation = async () => {
    if (!simulationId) return
    
    try {
      await api.pauseSimulation(simulationId)
      await updateStatus(simulationId)
      await refreshMetrics(simulationId)
      // accumulate runtime and stop clock
      if (runtimeStartMs !== null) {
        setRuntimeAccumulatedMs(prev => prev + (Date.now() - runtimeStartMs))
        setRuntimeStartMs(null)
      }
    } catch (error) {
      console.error('Failed to pause simulation:', error)
    }
  }

  const resumeSimulation = async () => {
    if (!simulationId) return
    
    try {
      await api.resumeSimulation(simulationId)
      await updateStatus(simulationId)
      // resume wall-clock
      if (runtimeStartMs === null) {
        setRuntimeStartMs(Date.now())
      }
    } catch (error) {
      console.error('Failed to resume simulation:', error)
    }
  }

  const stopSimulation = async () => {
    if (!simulationId) return
    
    try {
      await api.stopSimulation(simulationId)
      wsClient.disconnect()
      setIsConnected(false)
      await updateStatus(simulationId)
      // stop wall-clock without resetting accumulated (keeps total)
      if (runtimeStartMs !== null) {
        setRuntimeAccumulatedMs(prev => prev + (Date.now() - runtimeStartMs))
        setRuntimeStartMs(null)
      }
    } catch (error) {
      console.error('Failed to stop simulation:', error)
    }
  }

  const resetSimulation = async () => {
    if (!simulationId) return
    
    try {
      await api.resetSimulation(simulationId)
      await updateStatus(simulationId)
      // reset wall-clock
      setRuntimeAccumulatedMs(0)
      setRuntimeStartMs(Date.now())
    } catch (error) {
      console.error('Failed to reset simulation:', error)
    }
  }

  const updateStatus = async (id: string) => {
    try {
      const statusData = await api.getSimulationStatus(id)
      setStatus(statusData)
    } catch (error: any) {
      // If simulation not found (404), clear simulation ID to prevent repeated errors
      if (error?.message?.includes('404') || error?.message?.includes('not found')) {
        console.warn(`Simulation ${id} no longer exists. Clearing simulation ID.`)
        setSimulationId(null)
        setStatus(null)
      } else {
        console.error('Failed to update status:', error)
      }
    }
  }

  const refreshMetrics = async (id: string) => {
    try {
      const response = await api.getMetrics(id)
      setMetricsSnapshot(response.metrics)
    } catch (error) {
      console.error('Failed to refresh metrics:', error)
    }
  }

  const getStatusIndicator = () => {
    if (!status && !simulationData) return 'disconnected'
    const running = status?.is_running ?? true
    const paused = status?.is_paused ?? false
    if (running && !paused) return 'running'
    if (paused) return 'paused'
    return 'stopped'
  }

  const getStatusText = () => {
    if (!status && !simulationData) return 'Disconnected'
    const running = status?.is_running ?? true
    const paused = status?.is_paused ?? false
    if (running && !paused) return 'Running'
    if (paused) return 'Paused'
    return 'Stopped'
  }

  return (
    <div className="app">
      <header className="header">
        <div style={{ display: 'flex', alignItems: 'center', gap: '1rem' }}>
          <img src={logoImage} alt="Live 2.0" style={{ height: '160px', width: 'auto' }} />
          <div className="status">
            <div className={`status-indicator ${getStatusIndicator()}`} />
            <span>{getStatusText()}</span>
            {simulationId && (
              <span className="text-sm text-gray-400 ml-2">
                ID: {simulationId.slice(-8)}
              </span>
            )}
          </div>
        </div>
        
        <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', flex: 1 }}>
          <img src={headerImage} alt="Live 2.0 Header" style={{ height: '160px', width: 'auto' }} />
        </div>
        
        <div className="flex items-center gap-2">
          <button
            className={`btn ${view === 'api-v1' ? 'btn-primary' : 'btn-secondary'}`}
            onClick={() => setView(view === 'api-v1' ? 'simulation' : 'api-v1')}
            title="Toggle between Simulation and API v1 Jobs view"
          >
            <Database size={16} />
            {view === 'api-v1' ? 'Simulation' : 'API v1 Jobs'}
          </button>
          
          {view === 'simulation' && (
            <>
              <label htmlFor="mode-select" className="text-sm text-gray-300">Mode:</label>
              <select
                id="mode-select"
                value={mode}
                onChange={(e) => setMode(e.target.value as any)}
                className="text-sm px-2 py-1 border border-white/20 rounded bg-white/10 text-white"
                aria-label="Simulation mode selection"
              >
                <option value="open_chemistry">Open Chemistry</option>
                <option value="preset_prebiotic">Preset Prebiotic</option>
              </select>
              
              <button
                className="btn btn-primary"
                onClick={status?.is_running ? (status?.is_paused ? resumeSimulation : pauseSimulation) : (simulationId ? () => startSimulation(simulationId) : initializeSimulation)}
                disabled={false}
              >
                {status?.is_running ? (status?.is_paused ? <Play size={16} /> : <Pause size={16} />) : <Play size={16} />}
                {status?.is_running ? (status?.is_paused ? 'Resume' : 'Pause') : (simulationId ? 'Start' : 'Create & Start')}
              </button>
              
              <button
                className="btn btn-danger"
                onClick={stopSimulation}
                disabled={!simulationId}
              >
                <Square size={16} />
                Stop
              </button>
              
              <button
                className="btn btn-secondary"
                onClick={resetSimulation}
                disabled={!simulationId}
              >
                <RotateCcw size={16} />
                Reset
              </button>
            </>
          )}
        </div>
      </header>

      {view === 'api-v1' ? (
        <main className="main-content" style={{ padding: '2rem' }}>
          <APIv1Jobs apiBaseUrl="http://localhost:8001" />
        </main>
      ) : (
        <>
          <main className="main-content">
            <aside className="sidebar">
              <Controls
                simulationId={simulationId}
                status={status}
                onStatusUpdate={updateStatus}
                availableSpecies={availableSpecies}
                selectedSubstance={selectedSubstance}
                onSubstanceChange={setSelectedSubstance}
                runtimeMs={runtimeAccumulatedMs + (runtimeStartMs ? (runtimeNowMs - runtimeStartMs) : 0)}
              />
              
              <PerformancePanel
                simulationId={simulationId}
                performanceData={performanceData}
              />
              
              <GraphPreview
                data={simulationData}
                selectedSubstance={selectedSubstance}
                onSubstanceSelect={setSelectedSubstance}
              />
              
              {/* Novelty Panel with PubChem Matcher */}
              {simulationId && (
                <div className="mt-6">
                  <NoveltyPanel
                    simulationId={simulationId}
                    onSubstanceSelect={setSelectedSubstance}
                  />
                </div>
              )}
            </aside>

            <div className="canvas-container">
              <HeatmapCanvas
                data={simulationData}
                selectedSubstance={selectedSubstance}
                isConnected={isConnected}
              />
            </div>
          </main>
          
          {/* Legend at bottom */}
          <footer className="legend-footer" key={`legend-${simulationData?.step_count ?? status?.step_count ?? 0}`}>
            <div className="legend">
          <h3>🧪 Symulacja Molekularna Live 2.0</h3>
          <div className="legend-grid">
            <div className="legend-column">
              <h4>Elementy Wizualizacji</h4>
              <div className="legend-item">
                <div className="legend-dot" style={{ background: '#60a5fa' }}></div>
                <span>Cząstki molekularne</span>
              </div>
              <div className="legend-item">
                <div className="legend-dot" style={{ background: '#fbbf24' }}></div>
                <span>Pole energii (impulsy)</span>
              </div>
              <div className="legend-item">
                <div className="legend-dot" style={{ background: '#34d399' }}></div>
                <span>Wiązania chemiczne</span>
              </div>
              <div className="legend-item">
                <div className="legend-dot" style={{ background: '#f87171' }}></div>
                <span>Wysoka energia</span>
              </div>
            </div>
            <div className="legend-column" key={`status-${simulationData?.step_count ?? status?.step_count ?? 0}`}>
              <h4>Status Symulacji</h4>
              <div className="legend-info">
                <strong>Tryb:</strong> {mode === 'open_chemistry' ? 'Open Chemistry' : 'Preset Prebiotic'}
              </div>
              <div className="legend-info" key={`step-${simulationData?.step_count ?? status?.step_count ?? 0}`}>
                <strong>Krok:</strong> {simulationData?.step_count ?? status?.step_count ?? 0}
              </div>
              <div className="legend-info" key={`time-${simulationData?.current_time ?? 0}`}>
                <strong>Czas:</strong> {simulationData?.current_time?.toFixed(2) ?? '0.00'}s
              </div>
              <div className="legend-info">
                <strong>Cząstki:</strong> {metricsSnapshot?.particle_count ?? status?.particle_count ?? 0}
              </div>
            </div>
            <div className="legend-column">
              <h4>Co Obserwować</h4>
              <div className="legend-info">
                • <strong>Formowanie wiązań</strong> - cząstki łączą się w struktury
              </div>
              <div className="legend-info">
                • <strong>Impulsy energii</strong> - jasne plamy co 48 kroków
              </div>
              <div className="legend-info">
                • <strong>Samoorganizacja</strong> - spontaniczne tworzenie klastrów
              </div>
              <div className="legend-info">
                • <strong>Nowe substancje</strong> - wykrywanie w NoveltyPanel
              </div>
            </div>
          </div>
        </div>
        </footer>
        </>
      )}
    </div>
  )
}

export default App
