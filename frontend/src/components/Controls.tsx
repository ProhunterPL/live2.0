import React, { useState, useEffect } from 'react'
import { BarChart3, Download, Upload, Settings } from 'lucide-react'
import { SimulationAPI } from '../lib/api'
import type { SimulationStatus, Metrics } from '../lib/types'

interface ControlsProps {
  simulationId: string | null
  status: SimulationStatus | null
  onStatusUpdate: (id: string) => void
  availableSpecies?: string[]
  selectedSubstance?: string | null
  onSubstanceChange?: (id: string | null) => void
  runtimeMs?: number
  metricsOverride?: Metrics | null
  currentStep?: number | null
}

const Controls: React.FC<ControlsProps> = ({
  simulationId,
  status,
  onStatusUpdate,
  availableSpecies,
  selectedSubstance,
  onSubstanceChange,
  runtimeMs,
  metricsOverride,
  currentStep
}) => {
  const [metrics, setMetrics] = useState<Metrics | null>(null)
  const [novelSubstances, setNovelSubstances] = useState<any[]>([])
  const [showSettings, setShowSettings] = useState(false)
  const [snapshotFilename, setSnapshotFilename] = useState('')
  // species selection handled via props

  const api = new SimulationAPI()

  useEffect(() => {
    if (simulationId) {
      updateMetrics()
      updateNovelSubstances()
      // Poll concentrations species list if present via metrics hint soon
      
      // Update metrics every 5 seconds
      const interval = setInterval(() => {
        updateMetrics()
        updateNovelSubstances()
      }, 5000)

      return () => clearInterval(interval)
    }
  }, [simulationId])

  const updateMetrics = async () => {
    if (!simulationId) return
    
    try {
      const response = await api.getMetrics(simulationId)
      setMetrics(response.metrics)
    } catch (error: any) {
      // Stop polling if simulation doesn't exist
      if (error?.message?.includes('404') || error?.message?.includes('not found')) {
        console.warn('Simulation not found, stopping metrics updates')
        return
      }
      console.error('Failed to update metrics:', error)
    }
  }

  const updateNovelSubstances = async () => {
    if (!simulationId) return
    
    try {
      const response = await api.getNovelSubstances(simulationId, 5)
      setNovelSubstances(response.substances)
    } catch (error: any) {
      // Stop polling if simulation doesn't exist
      if (error?.message?.includes('404') || error?.message?.includes('not found')) {
        console.warn('Simulation not found, stopping novel substances updates')
        return
      }
      console.error('Failed to update novel substances:', error)
    }
  }

  const handleSaveSnapshot = async () => {
    if (!simulationId) return
    
    try {
      const filename = snapshotFilename || `snapshot_${Date.now()}.json`
      await api.saveSnapshot(simulationId, filename)
      alert('Snapshot saved successfully!')
      setSnapshotFilename('')
    } catch (error) {
      console.error('Failed to save snapshot:', error)
      alert('Failed to save snapshot')
    }
  }

  const handleLoadSnapshot = async () => {
    if (!simulationId || !snapshotFilename) return
    
    try {
      await api.loadSnapshot(simulationId, snapshotFilename)
      alert('Snapshot loaded successfully!')
      onStatusUpdate(simulationId)
    } catch (error) {
      console.error('Failed to load snapshot:', error)
      alert('Failed to load snapshot')
    }
  }

  const formatNumber = (num: number): string => {
    if (num >= 1000000) {
      return (num / 1000000).toFixed(1) + 'M'
    } else if (num >= 1000) {
      return (num / 1000).toFixed(1) + 'K'
    }
    return num.toFixed(1)
  }

  const formatTime = (seconds: number): string => {
    if (!Number.isFinite(seconds)) return '0:00.00'
    if (seconds < 60) {
      const mins = Math.floor(seconds / 60)
      const secs = Math.floor(seconds % 60)
      const ms = Math.floor((seconds - Math.floor(seconds)) * 100)
      return `${mins}:${secs.toString().padStart(2, '0')}.${ms.toString().padStart(2, '0')}`
    }
    const hours = Math.floor(seconds / 3600)
    const minutes = Math.floor((seconds % 3600) / 60)
    const secs = Math.floor(seconds % 60)
    return `${hours}:${minutes.toString().padStart(2, '0')}:${secs.toString().padStart(2, '0')}`
  }

  return (
    <div className="controls">
      <div className="control-group">
        <h2 className="text-lg font-semibold text-white mb-4">Simulation Controls</h2>
        
        {/* Status */}
        <div className="control-group">
          <label>Status</label>
          <div className="flex items-center gap-2">
            <div className={`w-2 h-2 rounded-full ${
              status?.is_running ? (status.is_paused ? 'bg-yellow-500' : 'bg-green-500') : 'bg-red-500'
            }`} />
            <span className="text-sm">
              {status?.is_running ? (status.is_paused ? 'Paused' : 'Running') : 'Stopped'}
            </span>
          </div>
        </div>

        {/* Time (wall-clock runtime) */}
        {typeof runtimeMs === 'number' && (
          <div className="control-group">
            <label>Simulation Time</label>
            <div className="text-sm font-mono">
              {formatTime((runtimeMs as number) / 1000)}
            </div>
          </div>
        )}

        {/* Step Count */}
        {status && (
          <div className="control-group">
            <label>Step Count</label>
            <div className="text-sm font-mono">
              {formatNumber(status.step_count)}
            </div>
          </div>
        )}
      </div>

      {/* Preset Mode: Species Selection (if concentrations present via frontend main view) */}
      {Array.isArray(availableSpecies) && (availableSpecies as string[]).length > 0 && (
        <div className="control-group">
          <label>Select Species (Preset Mode)</label>
          <select
            className="w-full text-sm"
            value={selectedSubstance || ''}
            onChange={(e) => onSubstanceChange?.(e.target.value || null)}
          >
            <option value="">All</option>
            {(availableSpecies as string[]).map((s: string) => (
              <option key={s} value={s}>{s}</option>
            ))}
          </select>
        </div>
      )}

      {/* Metrics */}
      {metrics && (
        <div className="control-group">
          <h3 className="text-md font-semibold text-white mb-2 flex items-center gap-2">
            <BarChart3 size={16} />
            Metrics
          </h3>
          <div className="metrics">
            <div className="metric">
              <div className="metric-label">Particles</div>
              <div className="metric-value">{formatNumber(metrics.particle_count ?? 0)}</div>
            </div>
            <div className="metric">
              <div className="metric-label">Bonds</div>
              <div className="metric-value">{formatNumber(metrics.bond_count ?? 0)}</div>
            </div>
            <div className="metric">
              <div className="metric-label">Clusters</div>
              <div className="metric-value">{formatNumber(metrics.cluster_count ?? 0)}</div>
            </div>
            <div className="metric">
              <div className="metric-label">Novelty Rate</div>
              <div className="metric-value">{(((metrics.novelty_rate ?? 0) * 100)).toFixed(1)}%</div>
            </div>
            <div className="metric">
              <div className="metric-label">Health Score</div>
              <div className="metric-value">{
                Number.isFinite(metrics.health_score ?? 0)
                  ? (((metrics.health_score ?? 0) * 100).toFixed(1) + '%')
                  : '0.0%'
              }</div>
            </div>
            <div className="metric">
              <div className="metric-label">Total Energy</div>
              <div className="metric-value">{formatNumber(metrics.total_energy ?? 0)}</div>
            </div>
          </div>
        </div>
      )}

      {/* Novel Substances */}
      {novelSubstances.length > 0 && (
        <div className="control-group">
          <h3 className="text-md font-semibold text-white mb-2">Recent Discoveries</h3>
          <div className="substance-list">
            {novelSubstances.slice(0, 3).map((substance) => (
              <div key={substance.id} className="substance-item">
                <div className="substance-id">{substance.id.slice(-8)}</div>
                <div className="substance-props">
                  Size: {substance.size} | Complexity: {substance.complexity.toFixed(1)}
                </div>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Snapshot Controls */}
      <div className="control-group">
        <h3 className="text-md font-semibold text-white mb-2 flex items-center gap-2">
          <Settings size={16} />
          Snapshots
        </h3>
        
        <div className="flex flex-col gap-2">
          <input
            type="text"
            placeholder="Snapshot filename"
            value={snapshotFilename}
            onChange={(e) => setSnapshotFilename(e.target.value)}
            className="text-sm"
          />
          
          <div className="flex gap-2">
            <button
              className="btn btn-secondary flex-1"
              onClick={handleSaveSnapshot}
              disabled={!simulationId}
            >
              <Download size={14} />
              Save
            </button>
            
            <button
              className="btn btn-secondary flex-1"
              onClick={handleLoadSnapshot}
              disabled={!simulationId || !snapshotFilename}
            >
              <Upload size={14} />
              Load
            </button>
          </div>
        </div>
      </div>

      {/* Settings Toggle */}
      <div className="control-group">
        <button
          className="btn btn-secondary w-full"
          onClick={() => setShowSettings(!showSettings)}
        >
          <Settings size={16} />
          {showSettings ? 'Hide' : 'Show'} Settings
        </button>
      </div>

      {/* Advanced Settings */}
      {showSettings && (
        <div className="control-group">
          <h3 className="text-md font-semibold text-white mb-2">Advanced Settings</h3>
          
          <div className="space-y-2">
            <div className="flex items-center justify-between">
              <label className="text-sm">Auto-save snapshots</label>
              <input type="checkbox" className="w-4 h-4" />
            </div>
            
            <div className="flex items-center justify-between">
              <label className="text-sm">Show particle IDs</label>
              <input type="checkbox" className="w-4 h-4" />
            </div>
            
            <div className="flex items-center justify-between">
              <label className="text-sm">High contrast mode</label>
              <input type="checkbox" className="w-4 h-4" />
            </div>
          </div>
        </div>
      )}
    </div>
  )
}

export default Controls
