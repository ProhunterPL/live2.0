import React, { useState, useEffect } from 'react'
import type { SimulationData } from '../lib/types'

interface GraphPreviewProps {
  data: SimulationData | null
  selectedSubstance: string | null
  onSubstanceSelect: (id: string | null) => void
}

const GraphPreview: React.FC<GraphPreviewProps> = ({
  data,
  selectedSubstance,
  onSubstanceSelect
}) => {
  const [previewData, setPreviewData] = useState<{
    particles: [number, number][]
    bonds: [number, number][]
    attributes: [number, number, number, number][]
  } | null>(null)

  useEffect(() => {
    if (data && data.particles && data.bonds) {
      // Create a preview of the largest cluster
      const clustersAny = data.clusters || []
      const clusters = (clustersAny as any[]).map((c) => Array.isArray(c) ? { particles: c } : c)
      if (clusters.length > 0) {
        // Find the largest cluster
        const largestCluster = clusters.reduce((max, cluster) => 
          cluster.particles.length > max.particles.length ? cluster : max
        )

        if (largestCluster.particles.length > 1) {
          // Get particles and bonds for this cluster
          const clusterParticles = largestCluster.particles
          const bondsAny = data.bonds as any[]
          const clusterBonds = bondsAny
            .map((b: any) => Array.isArray(b) ? { particle_i: b[0], particle_j: b[1] } : b)
            .filter((bond: any) => clusterParticles.includes(bond.particle_i) && clusterParticles.includes(bond.particle_j))

          // Map to local indices
          const particleMap = new Map<number, number>()
          ;(clusterParticles as number[]).forEach((particleId: number, index: number) => {
            particleMap.set(particleId, index)
          })

          const localBonds: [number, number][] = clusterBonds.map(bond => [
            particleMap.get(bond.particle_i) || 0,
            particleMap.get(bond.particle_j) || 0
          ])

          const localParticles: [number, number][] = (clusterParticles as number[]).map((particleId: number) => 
            (data.particles as any).positions[particleId] || [0, 0]
          )

          const localAttributes: [number, number, number, number][] = (clusterParticles as number[]).map((particleId: number) => 
            (data.particles as any).attributes[particleId] || [1, 0, 0, 0]
          )

          setPreviewData({
            particles: localParticles,
            bonds: localBonds,
            attributes: localAttributes
          })
        }
      }
    }
  }, [data])

  const renderGraph = () => {
    if (!previewData || previewData.particles.length === 0) {
      return (
        <div className="flex items-center justify-center h-full text-gray-400">
          No clusters to display
        </div>
      )
    }

    const { particles, bonds, attributes } = previewData
    const width = 200
    const height = 200
    const padding = 20

    // Calculate bounds
    const minX = Math.min(...particles.map(p => p[0]))
    const maxX = Math.max(...particles.map(p => p[0]))
    const minY = Math.min(...particles.map(p => p[1]))
    const maxY = Math.max(...particles.map(p => p[1]))

    const scaleX = (width - 2 * padding) / Math.max(maxX - minX, 1)
    const scaleY = (height - 2 * padding) / Math.max(maxY - minY, 1)
    const scale = Math.min(scaleX, scaleY)

    const centerX = width / 2
    const centerY = height / 2

    return (
      <svg width={width} height={height} className="graph-svg">
        {/* Background */}
        <rect width={width} height={height} fill="rgba(255, 255, 255, 0.05)" />
        
        {/* Bonds */}
        {bonds.map((bond, index) => {
          const [i, j] = bond
          const pos1 = particles[i]
          const pos2 = particles[j]
          
          if (!pos1 || !pos2) return null

          const x1 = centerX + (pos1[0] - (minX + maxX) / 2) * scale
          const y1 = centerY + (pos1[1] - (minY + maxY) / 2) * scale
          const x2 = centerX + (pos2[0] - (minX + maxX) / 2) * scale
          const y2 = centerY + (pos2[1] - (minY + maxY) / 2) * scale

          return (
            <line
              key={index}
              x1={x1}
              y1={y1}
              x2={x2}
              y2={y2}
              stroke="rgba(255, 255, 255, 0.4)"
              strokeWidth="1"
            />
          )
        })}

        {/* Particles */}
        {particles.map((pos, index) => {
          const attr = attributes[index]
          if (!attr) return null

          const [mass, chargeX, chargeY, chargeZ] = attr
          const x = centerX + (pos[0] - (minX + maxX) / 2) * scale
          const y = centerY + (pos[1] - (minY + maxY) / 2) * scale

          // Color based on charge
          const chargeMagnitude = Math.sqrt(chargeX * chargeX + chargeY * chargeY + chargeZ * chargeZ)
          const hue = chargeMagnitude > 0 ? (chargeX > 0 ? 0 : 120) : 240
          
          const radius = Math.max(2, mass * 3)

          return (
            <circle
              key={index}
              cx={x}
              cy={y}
              r={radius}
              fill={`hsl(${hue}, 70%, 60%)`}
              stroke={`hsl(${hue}, 70%, 80%)`}
              strokeWidth="0.5"
            />
          )
        })}
      </svg>
    )
  }

  const getClusterStats = () => {
    if (!previewData) return null

    const { particles, bonds, attributes } = previewData
    const totalMass = attributes.reduce((sum, attr) => sum + attr[0], 0)
    const avgMass = totalMass / attributes.length
    const density = bonds.length / (particles.length * (particles.length - 1) / 2)

    return {
      size: particles.length,
      bonds: bonds.length,
      density: density.toFixed(3),
      avgMass: avgMass.toFixed(2)
    }
  }

  const stats = getClusterStats()

  return (
    <div className="graph-preview">
      <h3 className="text-md font-semibold text-white mb-2">Largest Cluster</h3>
      
      {previewData && (
        <div className="mb-2">
          {renderGraph()}
          
          {stats && (
            <div className="mt-2 text-xs text-gray-300">
              <div>Size: {stats.size} particles</div>
              <div>Bonds: {stats.bonds}</div>
              <div>Density: {stats.density}</div>
              <div>Avg Mass: {stats.avgMass}</div>
            </div>
          )}
        </div>
      )}

      {/* Cluster Selection */}
      {data && data.clusters && data.clusters.length > 0 && (
        <div className="mt-4">
          <h4 className="text-sm font-semibold text-white mb-2">Select Cluster</h4>
          <div className="space-y-1">
            {(data.clusters as any[]).slice(0, 5).map((cluster: any, index: number) => (
              <button
                key={index}
                className={`w-full text-left px-2 py-1 rounded text-xs ${
                  selectedSubstance === `cluster_${index}` 
                    ? 'bg-blue-500 bg-opacity-20 text-blue-300' 
                    : 'bg-gray-500 bg-opacity-10 text-gray-300 hover:bg-gray-500 bg-opacity-20'
                }`}
                onClick={() => onSubstanceSelect(`cluster_${index}`)}
              >
                Cluster {index + 1}: {(Array.isArray(cluster) ? cluster.length : cluster.particles.length)} particles
              </button>
            ))}
          </div>
        </div>
      )}

      {/* Instructions */}
      <div className="mt-4 text-xs text-gray-400">
        <p>Hover over particles in the main view to see detailed information.</p>
        <p>Click on clusters to highlight them.</p>
      </div>
    </div>
  )
}

export default GraphPreview
