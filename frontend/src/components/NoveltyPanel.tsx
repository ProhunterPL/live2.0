import React, { useState, useEffect, useRef } from 'react'
import { TrendingUp, TestTube, Beaker, Microscope, Download } from 'lucide-react'
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
  onSave?: (id: string) => void
}

const NovelSubstanceCard: React.FC<NovelSubstanceProps> = ({ substance, onSelect, onSave }) => {
  const timeStr = new Date(substance.timestamp * 1000).toLocaleTimeString()
  
  const handleSaveClick = (e: React.MouseEvent) => {
    e.stopPropagation()
    onSave?.(substance.id)
  }
  
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
        <div className="flex items-center gap-2">
          <span className="text-xs text-gray-400">{timeStr}</span>
          <button
            onClick={handleSaveClick}
            className="p-1 hover:bg-gray-700 rounded transition-colors"
            title="Save cluster for PubChem matching"
          >
            <Download className="w-3 h-3 text-green-400" />
          </button>
        </div>
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
  const [isMatching, setIsMatching] = useState(false)

  const api = useRef(new SimulationAPI()).current

  useEffect(() => {
    if (!simulationId) {
      console.log('NoveltyPanel: No simulationId, skipping update')
      return
    }

    const updateNovelty = async () => {
      console.log('NoveltyPanel: Updating novelty data for simulation', simulationId)
      setLoading(true)
      try {
        // Get novel substances
        const novelResponse = await api.getNovelSubstances(simulationId, 20)
        console.log('NoveltyPanel: Received', novelResponse.substances?.length || 0, 'substances')
        // BUGFIX: Deduplicate substances by ID to prevent React key warnings
        const substances = novelResponse.substances || []
        const uniqueSubstances = substances.filter((substance: NovelSubstance, index: number, self: NovelSubstance[]) => 
          index === self.findIndex((s: NovelSubstance) => s.id === substance.id)
        )
        
        // FLICKER FIX: Only update if there are actual changes
        setNovelSubstances(prev => {
          // If lengths differ, update
          if (prev.length !== uniqueSubstances.length) {
            return uniqueSubstances
          }
          
          // Check if any IDs are different
          const prevIds = new Set(prev.map((s: NovelSubstance) => s.id))
          const newIds = new Set(uniqueSubstances.map((s: NovelSubstance) => s.id))
          const hasChanges = prev.some((s: NovelSubstance) => !newIds.has(s.id)) || 
                            uniqueSubstances.some((s: NovelSubstance) => !prevIds.has(s.id))
          
          // Only update if there are actual changes
          return hasChanges ? uniqueSubstances : prev
        })

        // Get metrics for novelty rate
        const metricsResponse = await api.getMetrics(simulationId)
        const metrics = metricsResponse.metrics
        
        setNoveltyRate(metrics.novelty_rate || 0)
        
        // Catalog stats are directly in metrics (not nested)
        setTotalNovel(metrics.total_novel || 0)
        setTotalDiscovered(metrics.total_discovered || 0)
      } catch (error: any) {
        // Stop polling if simulation doesn't exist
        if (error?.message?.includes('404') || error?.message?.includes('not found')) {
          console.warn('Simulation not found, stopping novelty updates')
          return
        }
        console.error('Failed to update novelty data:', error)
      } finally {
        setLoading(false)
      }
    }

    updateNovelty()
    
    // Update every 15 seconds (reduced frequency to minimize flicker)
    const interval = setInterval(updateNovelty, 15000)
    return () => clearInterval(interval)
  }, [simulationId, api])

  const handleMatchAllClusters = async () => {
    if (!simulationId || novelSubstances.length === 0) return
    
    setIsMatching(true)
    try {
      console.log(`🔍 Matching ${novelSubstances.length} clusters with PubChem...`)
      
      let successCount = 0
      for (const substance of novelSubstances) {
        try {
          await api.matchSubstanceToPubchem(simulationId, substance.id)
          successCount++
        } catch (error) {
          console.error(`Failed to match cluster ${substance.id}:`, error)
        }
      }
      
      alert(`✅ Matched ${successCount}/${novelSubstances.length} clusters with PubChem!\n\nCheck the matches/ directory for results.`)
    } catch (error) {
      console.error('Failed to match clusters:', error)
      alert('❌ Failed to match clusters with PubChem')
    } finally {
      setIsMatching(false)
    }
  }

  const handleSaveCluster = async (substanceId: string) => {
    if (!simulationId) return
    
    try {
      console.log(`🔍 Matching cluster ${substanceId} to PubChem...`)
      
      // Call backend API to match substance
      const result = await api.matchSubstanceToPubchem(simulationId, substanceId)
      
      if (result.success) {
        console.log('✅ Match completed!')
        console.log(`SMILES: ${result.smiles}`)
        
        if (result.pubchem_match) {
          console.log(`🎯 PubChem Match: CID ${result.pubchem_match.cid}`)
          console.log(`   Name: ${result.pubchem_match.name}`)
          console.log(`   Formula: ${result.pubchem_match.formula}`)
        } else {
          console.log('ℹ️  No PubChem match found')
        }
        
        console.log('\n📁 Generated files:')
        Object.entries(result.files).forEach(([key, path]) => {
          console.log(`   ${key}: ${path}`)
        })
        
        // Show success message to user
        const matchInfo = result.pubchem_match 
          ? `Matched to ${result.pubchem_match.name} (CID ${result.pubchem_match.cid})`
          : 'No PubChem match found'
        
        alert(`✅ Cluster matched successfully!\n\n${matchInfo}\n\nSMILES: ${result.smiles}\n\nFiles saved to matches/ folder.`)
      }
      
    } catch (error) {
      console.error('❌ Failed to match cluster:', error)
      alert(`Failed to match cluster: ${error}`)
    }
  }

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
    <div className="novelty-panel space-y-4 border-2 border-purple-500 rounded-lg p-4">
      {/* Header */}
      <div className="flex items-center gap-2 mb-4">
        <TrendingUp className="w-5 h-5 text-purple-400" />
        <h3 className="text-lg font-semibold text-white">🔬 PubChem Matcher</h3>
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
        <div className="flex items-center justify-between mb-3">
          <h4 className="text-md font-semibold text-white flex items-center gap-2">
            <TestTube className="w-4 h-4 text-green-400" />
            Recent Discoveries
          </h4>
          
          {/* Main PubChem Matcher Button */}
          {novelSubstances.length > 0 && (
            <button
              onClick={handleMatchAllClusters}
              disabled={isMatching}
              className="px-3 py-1 bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-700 hover:to-purple-700 text-white text-xs font-semibold rounded-lg shadow-lg transition-all disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2"
              title={`Match all ${novelSubstances.length} clusters with PubChem database`}
            >
              {isMatching ? (
                <>
                  <div className="animate-spin rounded-full h-3 w-3 border-b-2 border-white"></div>
                  <span>Matching...</span>
                </>
              ) : (
                <>
                  <Microscope className="w-3 h-3" />
                  <span>Match All ({novelSubstances.length})</span>
                </>
              )}
            </button>
          )}
        </div>
        
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
                onSave={handleSaveCluster}
              />
            ))}
          </div>
        ) : (
          <div className="text-center py-8 text-gray-400 bg-gray-800 bg-opacity-30 rounded border border-gray-600">
            <TestTube className="w-8 h-8 mx-auto mb-2 opacity-50" />
            <p className="font-semibold text-yellow-400">⏳ Waiting for clusters...</p>
            <p className="text-xs mt-2">No novel substances detected yet.</p>
            <p className="text-xs mt-1">Clusters will appear after 100+ simulation steps.</p>
            <p className="text-xs mt-1 text-blue-400">💡 Clusters need time to form through particle bonds.</p>
            <p className="text-xs mt-2 text-purple-400">Click Download 🔽 button to match with PubChem! 🧪</p>
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
