import React, { useState, useEffect } from 'react'
import { Activity, AlertTriangle, CheckCircle, Clock, Zap } from 'lucide-react'
import { SimulationAPI } from '../lib/api'

interface PerformanceMetrics {
  fps: number
  avg_step_time_ms: number
  min_step_time_ms: number
  max_step_time_ms: number
  avg_visualization_time_ms: number
  max_visualization_time_ms: number
  avg_broadcast_time_ms: number
  max_broadcast_time_ms: number
  performance_status: string
  performance_score: number
}

interface PerformancePanelProps {
  simulationId: string | null
  performanceData?: PerformanceMetrics | null
}

const PerformancePanel: React.FC<PerformancePanelProps> = ({
  simulationId,
  performanceData
}) => {
  const [performance, setPerformance] = useState<PerformanceMetrics | null>(null)
  const [warnings, setWarnings] = useState<string[]>([])
  const [lastUpdate, setLastUpdate] = useState<Date | null>(null)

  const api = new SimulationAPI()

  useEffect(() => {
    if (performanceData) {
      setPerformance(performanceData)
      setLastUpdate(new Date())
      return
    }

    if (!simulationId) return

    const updatePerformance = async () => {
      try {
        const response = await api.getPerformanceMetrics(simulationId)
        setPerformance(response.performance)
        setWarnings(response.warnings || [])
        setLastUpdate(new Date())
      } catch (error: any) {
        // Stop polling if simulation doesn't exist
        if (error?.message?.includes('404') || error?.message?.includes('not found')) {
          console.warn('Simulation not found, stopping performance metrics updates')
          return
        }
        console.error('Failed to update performance metrics:', error)
      }
    }

    updatePerformance()
    
    // Update every 5 seconds
    const interval = setInterval(updatePerformance, 5000)
    return () => clearInterval(interval)
  }, [simulationId, performanceData])

  const getStatusColor = (status: string) => {
    switch (status) {
      case 'excellent': return 'text-green-400'
      case 'good': return 'text-blue-400'
      case 'acceptable': return 'text-yellow-400'
      case 'poor': return 'text-red-400'
      case 'warming_up': return 'text-gray-400'
      default: return 'text-gray-400'
    }
  }

  const getStatusIcon = (status: string) => {
    switch (status) {
      case 'excellent': return <CheckCircle className="w-4 h-4" />
      case 'good': return <CheckCircle className="w-4 h-4" />
      case 'acceptable': return <AlertTriangle className="w-4 h-4" />
      case 'poor': return <AlertTriangle className="w-4 h-4" />
      case 'warming_up': return <Clock className="w-4 h-4" />
      default: return <Activity className="w-4 h-4" />
    }
  }

  const formatTime = (ms: number) => {
    if (ms < 1) return `${(ms * 1000).toFixed(1)}μs`
    if (ms < 1000) return `${ms.toFixed(1)}ms`
    return `${(ms / 1000).toFixed(2)}s`
  }

  if (!performance) {
    return (
      <div className="control-group">
        <h3 className="text-md font-semibold text-white mb-2 flex items-center gap-2">
          <Activity className="w-5 h-5 text-blue-400" />
          Performance Monitor
        </h3>
        <div className="text-sm text-gray-400">Loading performance data...</div>
      </div>
    )
  }

  return (
    <div className="control-group">
      <h3 className="text-md font-semibold text-white mb-2 flex items-center gap-2">
        <Activity className="w-5 h-5 text-blue-400" />
        Performance Monitor
        {lastUpdate && (
          <span className="text-xs text-gray-500 ml-auto">
            {lastUpdate.toLocaleTimeString()}
          </span>
        )}
      </h3>

      {/* Performance Status */}
      <div className="mb-3">
        <div className={`flex items-center gap-2 text-sm ${getStatusColor(performance.performance_status)}`}>
          {getStatusIcon(performance.performance_status)}
          <span className="capitalize">{performance.performance_status.replace('_', ' ')}</span>
          <span className="text-xs text-gray-400 ml-auto">
            Score: {(performance.performance_score * 100).toFixed(0)}%
          </span>
        </div>
      </div>

      {/* FPS */}
      <div className="mb-2">
        <div className="flex justify-between items-center text-sm">
          <span className="text-gray-300">FPS</span>
          <span className={`font-mono ${performance.fps >= 60 ? 'text-green-400' : performance.fps >= 30 ? 'text-yellow-400' : 'text-red-400'}`}>
            {performance.fps.toFixed(1)}
          </span>
        </div>
        <div className="w-full bg-gray-700 rounded-full h-1 mt-1">
          <div 
            className={`h-1 rounded-full ${performance.fps >= 60 ? 'bg-green-400' : performance.fps >= 30 ? 'bg-yellow-400' : 'bg-red-400'}`}
            style={{ width: `${Math.min(100, (performance.fps / 60) * 100)}%` }}
          ></div>
        </div>
      </div>

      {/* Step Times */}
      <div className="mb-2">
        <div className="text-xs text-gray-400 mb-1">Step Times</div>
        <div className="grid grid-cols-3 gap-2 text-xs">
          <div>
            <div className="text-gray-300">Avg</div>
            <div className="font-mono text-blue-400">{formatTime(performance.avg_step_time_ms)}</div>
          </div>
          <div>
            <div className="text-gray-300">Min</div>
            <div className="font-mono text-green-400">{formatTime(performance.min_step_time_ms)}</div>
          </div>
          <div>
            <div className="text-gray-300">Max</div>
            <div className="font-mono text-red-400">{formatTime(performance.max_step_time_ms)}</div>
          </div>
        </div>
      </div>

      {/* Visualization Times */}
      <div className="mb-2">
        <div className="text-xs text-gray-400 mb-1">Visualization</div>
        <div className="grid grid-cols-2 gap-2 text-xs">
          <div>
            <div className="text-gray-300">Avg</div>
            <div className="font-mono text-blue-400">{formatTime(performance.avg_visualization_time_ms)}</div>
          </div>
          <div>
            <div className="text-gray-300">Max</div>
            <div className="font-mono text-red-400">{formatTime(performance.max_visualization_time_ms)}</div>
          </div>
        </div>
      </div>

      {/* Broadcast Times */}
      <div className="mb-2">
        <div className="text-xs text-gray-400 mb-1">Broadcast</div>
        <div className="grid grid-cols-2 gap-2 text-xs">
          <div>
            <div className="text-gray-300">Avg</div>
            <div className="font-mono text-blue-400">{formatTime(performance.avg_broadcast_time_ms)}</div>
          </div>
          <div>
            <div className="text-gray-300">Max</div>
            <div className="font-mono text-red-400">{formatTime(performance.max_broadcast_time_ms)}</div>
          </div>
        </div>
      </div>

      {/* Warnings */}
      {warnings.length > 0 && (
        <div className="mt-3 p-2 bg-red-900 bg-opacity-20 border border-red-500 border-opacity-30 rounded">
          <div className="flex items-center gap-2 text-red-400 text-sm mb-1">
            <AlertTriangle className="w-4 h-4" />
            <span className="font-semibold">Performance Warnings</span>
          </div>
          {warnings.map((warning, index) => (
            <div key={index} className="text-xs text-red-300">
              • {warning}
            </div>
          ))}
        </div>
      )}

      {/* Performance Tips */}
      {performance.performance_status === 'poor' && (
        <div className="mt-3 p-2 bg-yellow-900 bg-opacity-20 border border-yellow-500 border-opacity-30 rounded">
          <div className="flex items-center gap-2 text-yellow-400 text-sm mb-1">
            <Zap className="w-4 h-4" />
            <span className="font-semibold">Performance Tips</span>
          </div>
          <div className="text-xs text-yellow-300">
            • Reduce particle count<br/>
            • Lower visualization frequency<br/>
            • Check system resources<br/>
            • Restart simulation if needed
          </div>
        </div>
      )}
    </div>
  )
}

export default PerformancePanel
