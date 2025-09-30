import React, { useState, useEffect, useRef } from 'react'
import { Play, Pause, Square, RotateCcw } from 'lucide-react'
import HeatmapCanvas from './components/HeatmapCanvas'
import Controls from './components/Controls'
import GraphPreview from './components/GraphPreview'
import { SimulationAPI } from './lib/api.ts'
import { WebSocketClient } from './lib/ws.ts'
import type { SimulationData, SimulationStatus } from './lib/types'

const App: React.FC = () => {
  const [simulationId, setSimulationId] = useState<string | null>(null)
  const [status, setStatus] = useState<SimulationStatus | null>(null)
  const [simulationData, setSimulationData] = useState<SimulationData | null>(null)
  const [isConnected, setIsConnected] = useState(false)
  const [selectedSubstance, setSelectedSubstance] = useState<string | null>(null)
  const [mode, setMode] = useState<'preset_prebiotic' | 'open_chemistry'>('open_chemistry')
  const [runtimeStartMs, setRuntimeStartMs] = useState<number | null>(null)
  const [runtimeAccumulatedMs, setRuntimeAccumulatedMs] = useState<number>(0)
  const [runtimeNowMs, setRuntimeNowMs] = useState<number>(Date.now())
  
  const api = useRef(new SimulationAPI()).current
  const wsClient = useRef(new WebSocketClient()).current

  const initializedRef = useRef(false)
  useEffect(() => {
    // Prevent duplicate init in React 18 StrictMode
    if (initializedRef.current) return
    initializedRef.current = true
    initializeSimulation()
    
    return () => {
      // Cleanup: only disconnect WS; don't stop simulation implicitly in dev
      wsClient.disconnect()
    }
  }, [wsClient])

  // Poll status periodically so the time advances even if WS data is delayed
  useEffect(() => {
    if (!simulationId) return
    const interval = setInterval(() => {
      updateStatus(simulationId)
    }, 1000)
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
          config: {
            grid_height: 128,
            grid_width: 128,
            mode,
            max_particles: 5000,
            max_time: 1000,
            dt: 0.01,
            energy_decay: 0.95,
            energy_threshold: 0.1,
            particle_radius: 0.5,
            binding_threshold: 0.8,
            unbinding_threshold: 0.2,
            novelty_window: 100,
            min_cluster_size: 2,
            vis_frequency: 5,
            log_frequency: 100,
            seed: Math.floor(Math.random() * 1000000)
          },
          mode
        })
        if (response.success && response.simulation_id) {
          setSimulationId(response.simulation_id)
          await connectWebSocket(response.simulation_id)
          // Don't auto-start simulation - wait for user to press play
          // await startSimulation(response.simulation_id)
          setRuntimeStartMs(null)
        }
      } catch (err) {
        console.error('Failed to recreate simulation:', err)
      }
    }
    recreate()
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [mode, api, wsClient])

  const initializeSimulation = async () => {
    try {
      const response = await api.createSimulation({
        config: {
          grid_height: 128,
          grid_width: 128,
          mode,
          max_particles: 5000,
          max_time: 1000,
          dt: 0.01,
          energy_decay: 0.99,
          energy_threshold: 0.01,
          particle_radius: 0.5,
          binding_threshold: 0.8,
          unbinding_threshold: 0.2,
          novelty_window: 100,
          min_cluster_size: 2,
            vis_frequency: 5,
          log_frequency: 100,
          seed: Math.floor(Math.random() * 1000000)
        },
        mode
      })

      if (response.success && response.simulation_id) {
        setSimulationId(response.simulation_id)
        await connectWebSocket(response.simulation_id)
        // Don't auto-start simulation - wait for user to press play
        // await startSimulation(response.simulation_id)
        // start wall-clock runtime
        setRuntimeAccumulatedMs(0)
        setRuntimeStartMs(null)
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
      // Pass only the simulationId; ws client will build the full URL
      await wsClient.connect(id)
      setIsConnected(true)
      
      wsClient.on('simulation_data', (data: SimulationData) => {
        setSimulationData(data)
      })
      wsClient.on('disconnected', () => {
        setIsConnected(false)
      })
      wsClient.on('error', (error: Error) => {
        console.error('WebSocket error:', error)
        setIsConnected(false)
      })
    } catch (error) {
      console.error('Failed to connect WebSocket:', error)
    }
  }

  const startSimulation = async (id: string) => {
    try {
      await api.startSimulation(id)
      await updateStatus(id)
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
    } catch (error) {
      console.error('Failed to update status:', error)
    }
  }

  const getStatusIndicator = () => {
    if (!status) return 'disconnected'
    if (status.is_running && !status.is_paused) return 'running'
    if (status.is_paused) return 'paused'
    return 'stopped'
  }

  const getStatusText = () => {
    if (!status) return 'Disconnected'
    if (status.is_running && !status.is_paused) return 'Running'
    if (status.is_paused) return 'Paused'
    return 'Stopped'
  }

  return (
    <div className="app">
      <header className="header">
        <div className="status">
          <div className={`status-indicator ${getStatusIndicator()}`} />
          <span>{getStatusText()}</span>
          {simulationId && (
            <span className="text-sm text-gray-400 ml-2">
              ID: {simulationId.slice(-8)}
            </span>
          )}
        </div>
        
        <div className="flex items-center gap-2">
          <label className="text-sm text-gray-300">Mode:</label>
          <select
            value={mode}
            onChange={(e) => setMode(e.target.value as any)}
            className="text-sm px-2 py-1 border border-white/20 rounded bg-white/10 text-white"
          >
            <option value="open_chemistry">Open Chemistry</option>
            <option value="preset_prebiotic">Preset Prebiotic</option>
          </select>
          
          <button
            className="btn btn-primary"
            onClick={status?.is_running ? (status?.is_paused ? resumeSimulation : pauseSimulation) : () => startSimulation(simulationId!)}
            disabled={!simulationId}
          >
            {status?.is_running ? (status?.is_paused ? <Play size={16} /> : <Pause size={16} />) : <Play size={16} />}
            {status?.is_running ? (status?.is_paused ? 'Resume' : 'Pause') : 'Start'}
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
        </div>
      </header>

      <main className="main-content">
        <aside className="sidebar">
          <Controls
            simulationId={simulationId}
            status={status}
            onStatusUpdate={updateStatus}
            runtimeMs={runtimeAccumulatedMs + (runtimeStartMs ? (runtimeNowMs - runtimeStartMs) : 0)}
          />
          
          {simulationData && (
            <GraphPreview
              data={simulationData}
              selectedSubstance={selectedSubstance}
              onSubstanceSelect={setSelectedSubstance}
            />
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
    </div>
  )
}

export default App
