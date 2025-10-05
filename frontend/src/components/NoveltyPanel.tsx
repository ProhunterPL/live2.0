import React, { useState, useEffect } from 'react'
import { TrendingUp, TestTube, Beaker, Microscope } from 'lucide-react'
import { SimulationAPI } from '../lib/api'

interface NovelSubstance {
  id: string
  timestamp: number  
  size: number
  complexity: number
  properties: {
    mass: number
    charge: [number, number, number]
    bonds: number
    graph_density: number
  }
}

interface NovelSubstanceProps {
  substance: NovelSubstance
  onSelect?: (id: string) => void
}

const NovelSubstanceCard: React.FC<NovelSubstanceProps> = ({ substance, onSelect }) => {
  const timeStr = new Date(substance.timestamp * 1000).toLocaleTimeString()
  
  return (
    <div 
      className="bg-gray-800 bg-opacity-50 p-3 rounded border border-gray-600 hover:border-blue-400 cursor-pointer transition-colors"
      onClick={() => onSelect?.(substance.id)}
    >
      <div className="flex items-center justify-between mb-2">
        <div className="flex items-center gap-2">
          <Beaker className="w-4 h-4 text-blue-400" />
          <span className="text-sm font-medium text-blue-300">{substance.id.slice(-8)}</span>
        </div>
        <span className="text-xs text-gray-400">{timeStr}</span>
      </div>
      
      <div className="space-y-1 text-xs text-gray-300">
        <div className="flex justify-between">
          <span>Size:</span>
          <span className="text-blue-300">{substance.size} atoms</span>
        </div>
        <div className="flex justify-between">
          <span>Complexity:</span>
          <span className="text-green-300">{substance.complexity.toFixed(3)}</span>
        </div>
        <div className="flex justify-between">
          <span>Bonds:</span>
          <span className="text-yellow-300">{substance.properties.bonds}</span>
        </div>
        <div className="flex justify-between">
          <span>Density:</span>
          <span className="text-purple-300">{substance.properties.graph_density.toFixed(3)}</span>
        </div>
        <div className="flex justify-between">
          <span>Avg Mass:</span>
          <span className="text-orange-300">{substance.properties.mass.toFixed(2)}</span>
        </div>
      </div>
    </div>
  )
}

interface NoveltyPanelProps {
  simulationId: string | null
  onSubstanceSelect?: (id: string) => void
}

const NoveltyPanel: React.FC<NoveltyPanelProps> = ({ simulationId, onSubstanceSelect }) => {
  const [novelSubstances, setNovelSubstances] = useState<NovelSubstance[]>([])
  const [noveltyRate, setNoveltyRate] = useState(0)
  const [totalNovel, setTotalNovel] = useState(0)
  const [totalDiscovered, setTotalDiscovered] = useState(0)
  const [loading, setLoading] = useState(false)

  const api = new SimulationAPI()

  useEffect(() => {
    if (!simulationId) return

    const updateNovelty = async () => {
      setLoading(true)
      try {
        // Get novel substances
        const novelResponse = await api.getNovelSubstances(simulationId, 20)
        setNovelSubstances(novelResponse.substances || [])

        // Get metrics for novelty rate
        const metricsResponse = await api.getMetrics(simulationId)
        const metrics = metricsResponse.metrics
        
        setNoveltyRate(metrics.novelty_rate || 0)
        
        // Parse catalog stats if available
        if (metrics.catalog_stats) {
          setTotalNovel(metrics.catalog_stats.total_novel || 0)
          setTotalDiscovered(metrics.catalog_stats.total_discovered || 0)
        }
      } catch (error) {
        console.error('Failed to update novelty data:', error)
      } finally {
        setLoading(false)
      }
    }

    updateNovelty()
    
    // Update every 10 seconds
    const interval = setInterval(updateNovelty, 10000)
    return () => clearInterval(interval)
  }, [simulationId, api])

  const getNoveltyColor = (rate: number) => {
    if (rate > 0.3) return 'text-red-400'
    if (rate > 0.1) return 'text-yellow-400'
    return 'text-green-400'
  }

  const getNoveltyMessage = (rate: number) => {
    if (rate > 0.3) return '🔥 High Emergence'
    if (rate > 0.1) return '⚡ Medium Activity'
    if (rate > 0.01) return '🌱 Low Activity'
    return '💤 Dormant'
  }

  return (
    <div className="novelty-panel space-y-4">
      {/* Header */}
      <div className="flex items-center gap-2 mb-4">
        <TrendingUp className="w-5 h-5 text-purple-400" />
        <h3 className="text-lg font-semibold text-white">Novelty Detection</h3>
      </div>

      {/* Novelty Stats */}
      <div className="bg-gray-800 bg-opacity-40 p-4 rounded border border-purple-500 border-opacity-30">
        <div className="space-y-3">
          <div className="flex items-center justify-between">
            <span className="text-sm text-gray-300">Emergence Rate:</span>
            <span className={`text-lg font-bold ${getNoveltyColor(noveltyRate)}`}>
              {" "}({noveltyRate.toFixed(4)})
            </span>
          </div>
          
          <div className="flex items-center gap-2">
            <Microscope className="w-4 h-4 text-blue-400" />
            <span className="text-sm text-gray-300">{getNoveltyMessage(noveltyRate)}</span>
          </div>
          
          <div className="grid grid-cols-2 gap-4 text-xs">
            <div className="text-center">
              <div className="text-blue-300 font-semibold">{totalNovel}</div>
              <div className="text-gray-400">Novel Substances</div>
            </div>
            <div className="text-center">
              <div className="text-green-300 font-semibold">{totalDiscovered}</div>
              <div className="text-gray-400">Total Discovered</div>
            </div>
          </div>
        </div>
      </div>

      {/* Recent Discoveries */}
      <div>
        <h4 className="text-md font-semibold text-white mb-3 flex items-center gap-2">
          <TestTube className="w-4 h-4 text-green-400" />
          Recent Discoveries
        </h4>
        
        {loading ? (
          <div className="flex items-center justify-center py-8">
            <div className="animate-spin rounded-full h-6 w-6 border-b-2 border-blue-400"></div>
          </div>
        ) : novelSubstances.length > 0 ? (
          <div className="space-y-2 max-h-64 overflow-y-auto">
            {novelSubstances.map((substance) => (
              <NovelSubstanceCard 
                key={substance.id}
                substance={substance}
                onSelect={onSubstanceSelect}
              />
            ))}
          </div>
        ) : (
          <div className="text-center py-8 text-gray-400">
            <TestTube className="w-8 h-8 mx-auto mb-2 opacity-50" />
            <p>No novel substances detected yet...</p>
            <p className="text-xs mt-1">Start simulation to detect emergence</p>
          </div>
        )}
      </div>

      {/* Instructions */}
      <div className="text-xs text-gray-400 bg-gray-800 bg-opacity-30 p-3 rounded">
        <p className="font-medium text-gray-300 mb-1">About Novelty:</p>
        <p>
          Novel substances are unique molecular structures discovered during simulation. 
          Higher emergence rates indicate active chemical evolution and new structure formation.
        </p>
      </div>
    </div>
  )
}

export default NoveltyPanel
