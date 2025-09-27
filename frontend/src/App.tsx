import React, { useState, useEffect } from 'react'
import { Play, Pause, Square, RotateCcw, Settings, BarChart3 } from 'lucide-react'
import HeatmapCanvas from './components/HeatmapCanvas'
import Controls from './components/Controls'
import GraphPreview from './components/GraphPreview'
import { SimulationAPI } from './lib/api'
import { WebSocketClient } from './lib/ws'
import type { SimulationData, SimulationStatus } from './lib/types'

const App: React.FC = () => {
  const [simulationId, setSimulationId] = useState<string | null>(null)
  const [status, setStatus] = useState<SimulationStatus | null>(null)
  const [simulationData, setSimulationData] = useState<SimulationData | null>(null)
  const [isConnected, setIsConnected] = useState(false)
  const [selectedSubstance, setSelectedSubstance] = useState<string | null>(null)
  
  const api = new SimulationAPI()
  const wsClient = new WebSocketClient()

  useEffect(() => {
    // Initialize simulation on mount
    initializeSimulation()
    
    return () => {
      // Cleanup on unmount
      if (simulationId) {
        api.stopSimulation(simulationId)
        wsClient.disconnect()
      }
    }
  }, [])

  const initializeSimulation = async () => {
    try {
      const response = await api.createSimulation({
        config: {
          grid_height: 256,
          grid_width: 256,
          mode: 'open_chemistry',
          max_particles: 10000,
          dt: 0.01,
          energy_decay: 0.95,
          energy_threshold: 0.1,
          particle_radius: 0.5,
          binding_threshold: 0.8,
          unbinding_threshold: 0.2,
          novelty_window: 100,
          min_cluster_size: 2,
          vis_frequency: 10,
          log_frequency: 100,
          seed: Math.floor(Math.random() * 1000000)
        },
        mode: 'open_chemistry'
      })

      if (response.success && response.simulation_id) {
        setSimulationId(response.simulation_id)
        await connectWebSocket(response.simulation_id)
        await startSimulation(response.simulation_id)
      }
    } catch (error) {
      console.error('Failed to initialize simulation:', error)
    }
  }

  const connectWebSocket = async (id: string) => {
    try {
      await wsClient.connect(`ws://localhost:8000/simulation/${id}/stream`)
      setIsConnected(true)
      
      wsClient.onMessage = (data: SimulationData) => {
        setSimulationData(data)
      }
      
      wsClient.onClose = () => {
        setIsConnected(false)
      }
      
      wsClient.onError = (error: Error) => {
        console.error('WebSocket error:', error)
        setIsConnected(false)
      }
    } catch (error) {
      console.error('Failed to connect WebSocket:', error)
    }
  }

  const startSimulation = async (id: string) => {
    try {
      await api.startSimulation(id)
      await updateStatus(id)
    } catch (error) {
      console.error('Failed to start simulation:', error)
    }
  }

  const pauseSimulation = async () => {
    if (!simulationId) return
    
    try {
      await api.pauseSimulation(simulationId)
      await updateStatus(simulationId)
    } catch (error) {
      console.error('Failed to pause simulation:', error)
    }
  }

  const resumeSimulation = async () => {
    if (!simulationId) return
    
    try {
      await api.resumeSimulation(simulationId)
      await updateStatus(simulationId)
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
    } catch (error) {
      console.error('Failed to stop simulation:', error)
    }
  }

  const resetSimulation = async () => {
    if (!simulationId) return
    
    try {
      await api.resetSimulation(simulationId)
      await updateStatus(simulationId)
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
          <button
            className="btn btn-primary"
            onClick={status?.is_paused ? resumeSimulation : pauseSimulation}
            disabled={!status?.is_running}
          >
            {status?.is_paused ? <Play size={16} /> : <Pause size={16} />}
            {status?.is_paused ? 'Resume' : 'Pause'}
          </button>
          
          <button
            className="btn btn-danger"
            onClick={stopSimulation}
            disabled={!status?.is_running}
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
