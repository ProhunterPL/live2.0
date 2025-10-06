import React, { useState, useEffect, useRef } from 'react'
import { Download } from 'lucide-react'
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
    energies: number[]
  } | null>(null)
  const svgRef = useRef<SVGSVGElement>(null)

  const exportAsImage = async (format: 'png' | 'jpg' = 'png') => {
    if (!svgRef.current || !previewData) return

    try {
      const svg = svgRef.current
      const svgData = new XMLSerializer().serializeToString(svg)
      const canvas = document.createElement('canvas')
      const ctx = canvas.getContext('2d')
      
      if (!ctx) return

      // Set canvas size (higher resolution for better quality)
      const scale = 4 // 4x resolution for 800x800px final image
      canvas.width = 200 * scale
      canvas.height = 200 * scale
      
      // Create image from SVG
      const img = new Image()
      const svgBlob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' })
      const url = URL.createObjectURL(svgBlob)
      
      img.onload = () => {
        // Fill background with white
        ctx.fillStyle = '#ffffff'
        ctx.fillRect(0, 0, canvas.width, canvas.height)
        
        // Draw SVG
        ctx.drawImage(img, 0, 0, canvas.width, canvas.height)
        
        // Add timestamp and simulation info
        const now = new Date()
        const timestamp = now.toLocaleString('pl-PL', {
          year: 'numeric',
          month: '2-digit',
          day: '2-digit',
          hour: '2-digit',
          minute: '2-digit',
          second: '2-digit'
        })
        
        // Add simulation info if available
        const simulationInfo = data ? `Live2.0 Simulation - ${data.step_count || 0} steps` : 'Live2.0 Simulation'
        
        // Style for text (fixed size without scaling)
        ctx.font = `10px Arial`
        ctx.textAlign = 'left'
        ctx.textBaseline = 'top'
        
        // Calculate text dimensions
        const timestampMetrics = ctx.measureText(timestamp)
        const infoMetrics = ctx.measureText(simulationInfo)
        const padding = 4
        const lineHeight = 12
        
        // Calculate background size
        const maxWidth = Math.max(timestampMetrics.width, infoMetrics.width)
        const bgWidth = maxWidth + padding * 2
        const bgHeight = lineHeight * 2 + padding * 2
        
        // Add background for text
        ctx.fillStyle = 'rgba(255, 255, 255, 0.9)'
        ctx.fillRect(15 * scale, canvas.height - bgHeight - 15 * scale, bgWidth, bgHeight)
        
        // Add timestamp text
        ctx.fillStyle = '#333333'
        ctx.fillText(timestamp, 15 * scale + padding, canvas.height - bgHeight - 15 * scale + padding)
        
        // Add simulation info
        ctx.fillStyle = '#666666'
        ctx.fillText(simulationInfo, 15 * scale + padding, canvas.height - bgHeight - 15 * scale + padding + lineHeight)
        
        // Add cluster parameters in top left corner
        if (previewData) {
          const stats = getClusterStats()
          if (stats) {
            const clusterInfo = [
              `Size: ${stats.size} particles`,
              `Bonds: ${stats.bonds}`,
              `Density: ${stats.density}`,
              `Avg Mass: ${stats.avgMass}`,
              `Total Energy: ${stats.totalEnergy}`,
              `Avg Energy: ${stats.avgEnergy}`
            ]
            
            // Calculate dimensions for cluster info
            const clusterMetrics = clusterInfo.map(text => ctx.measureText(text))
            const maxClusterWidth = Math.max(...clusterMetrics.map(m => m.width))
            const clusterPadding = 4
            const clusterLineHeight = 12
            const clusterBgWidth = maxClusterWidth + clusterPadding * 2
            const clusterBgHeight = clusterInfo.length * clusterLineHeight + clusterPadding * 2
            
            // Add background for cluster info
            ctx.fillStyle = 'rgba(255, 255, 255, 0.9)'
            ctx.fillRect(15 * scale, 15 * scale, clusterBgWidth, clusterBgHeight)
            
            // Add cluster info text
            ctx.fillStyle = '#333333'
            clusterInfo.forEach((line, index) => {
              ctx.fillText(line, 15 * scale + clusterPadding, 15 * scale + clusterPadding + index * clusterLineHeight)
            })
          }
        }
        
        // Add Live2.0 watermark
        const watermarkImg = new Image()
        watermarkImg.onload = () => {
          // Calculate watermark size (fixed 96x150px without scaling)
          const watermarkWidth = 96
          const watermarkHeight = 150
          const watermarkX = canvas.width - watermarkWidth - 15 * scale
          const watermarkY = 15 * scale
          
          // Add shadow for watermark
          ctx.shadowColor = 'rgba(0, 0, 0, 0.3)'
          ctx.shadowBlur = 4 * scale
          ctx.shadowOffsetX = 2 * scale
          ctx.shadowOffsetY = 2 * scale
          
          // Draw watermark with transparency
          ctx.globalAlpha = 0.8
          ctx.drawImage(watermarkImg, watermarkX, watermarkY, watermarkWidth, watermarkHeight)
          
          // Reset shadow and alpha
          ctx.shadowColor = 'transparent'
          ctx.shadowBlur = 0
          ctx.shadowOffsetX = 0
          ctx.shadowOffsetY = 0
          ctx.globalAlpha = 1.0
          
          // Convert to desired format and download
          const mimeType = format === 'png' ? 'image/png' : 'image/jpeg'
          const quality = format === 'jpg' ? 0.9 : undefined
          
          canvas.toBlob((blob) => {
            if (blob) {
              const link = document.createElement('a')
              link.download = `largest_cluster_${Date.now()}.${format}`
              link.href = URL.createObjectURL(blob)
              link.click()
              URL.revokeObjectURL(link.href)
            }
          }, mimeType, quality)
          
          URL.revokeObjectURL(url)
        }
        
        // Load watermark image
        watermarkImg.src = '/images/logolivebw.png'
      }
      
      img.src = url
    } catch (error) {
      console.error('Error exporting image:', error)
    }
  }

  useEffect(() => {
    if (data && data.particles && data.bonds) {
      const clustersAny = data.clusters || []
      const clusters = (clustersAny as any[]).map((c) => Array.isArray(c) ? { particles: c } : c)
      
      if (clusters.length > 0) {
        let targetCluster = null
        
        // Check if a specific cluster is selected
        if (selectedSubstance && selectedSubstance.startsWith('cluster_')) {
          const clusterIndex = parseInt(selectedSubstance.replace('cluster_', ''))
          if (clusterIndex >= 0 && clusterIndex < clusters.length) {
            targetCluster = clusters[clusterIndex]
          }
        }
        
        // Fallback to largest cluster if no specific selection or invalid index
        if (!targetCluster) {
          targetCluster = clusters.reduce((max, cluster) => 
            cluster.particles.length > max.particles.length ? cluster : max
          )
        }

        if (targetCluster.particles.length > 1) {
          // Get particles and bonds for this cluster
          const clusterParticles = targetCluster.particles
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

          const localEnergies: number[] = (clusterParticles as number[]).map((particleId: number) => 
            (data.particles as any).energies?.[particleId] || 0
          )

          setPreviewData({
            particles: localParticles,
            bonds: localBonds,
            attributes: localAttributes,
            energies: localEnergies
          })
        }
      }
    }
  }, [data, selectedSubstance])

  const renderGraph = () => {
    if (!previewData || previewData.particles.length === 0) {
      return (
        <div className="flex items-center justify-center h-full text-gray-400">
          No clusters to display
        </div>
      )
    }

    const { particles, bonds,attributes, energies } = previewData
    const width = 200
    const height = 200
    const padding = 30 // Increased padding to avoid overlap with text

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
      <svg ref={svgRef} width={width} height={height} className="graph-svg">
        {/* Background */}
        <rect width={width} height={height} fill="#1a1a1a" />
        
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
              stroke="#ffffff"
              strokeWidth="2"
              opacity="0.8"
            />
          )
        })}

        {/* Particles */}
        {particles.map((pos, index) => {
          const attr = attributes[index]
          const energy = energies[index] || 0
          if (!attr) return null

          const [mass] = attr
          const x = centerX + (pos[0] - (minX + maxX) / 2) * scale
          const y = centerY + (pos[1] - (minY + maxY) / 2) * scale

          // Color based on energy (brighter = more energy)
          const energyNormalized = Math.min(energy / 10, 1) // Normalize energy to 0-1
          const hue = 200 + energyNormalized * 120 // Blue to red based on energy
          const saturation = 60 + energyNormalized * 40 // More saturated for higher energy
          const lightness = 50 + energyNormalized * 30 // Brighter for higher energy
          
          const radius = Math.max(3, mass * 2 + energyNormalized * 2)

          return (
            <g key={index}>
              <circle
                cx={x}
                cy={y}
                r={radius}
                fill={`hsl(${hue}, ${saturation}%, ${lightness}%)`}
                stroke="#ffffff"
                strokeWidth="1"
              />
              {/* Energy text */}
              <text
                x={x}
                y={y + 1}
                textAnchor="middle"
                fontSize="8"
                fill="#ffffff"
                fontWeight="bold"
              >
                {energy.toFixed(1)}
              </text>
              {/* Mass text */}
              <text
                x={x}
                y={y + 10}
                textAnchor="middle"
                fontSize="6"
                fill="#cccccc"
              >
                {mass.toFixed(1)}
              </text>
            </g>
          )
        })}
      </svg>
    )
  }

  const getClusterStats = () => {
    if (!previewData) return null

    const { particles, bonds,attributes, energies } = previewData
    const totalMass = attributes.reduce((sum, attr) => sum + attr[0], 0)
    const avgMass = totalMass / attributes.length
    const totalEnergy = energies.reduce((sum, energy) => sum + energy, 0)
    const avgEnergy = totalEnergy / energies.length
    const density = bonds.length / (particles.length * (particles.length - 1) / 2)

    return {
      size: particles.length,
      bonds: bonds.length,
      density: density.toFixed(3),
      avgMass: avgMass.toFixed(2),
      totalEnergy: totalEnergy.toFixed(2),
      avgEnergy: avgEnergy.toFixed(2)
    }
  }

  const stats = getClusterStats()

  // Get current cluster info for title
  const getCurrentClusterInfo = () => {
    if (!data || !data.clusters) return null
    
    const clusters = (data.clusters as any[]).map((c) => Array.isArray(c) ? { particles: c } : c)
    
    if (selectedSubstance && selectedSubstance.startsWith('cluster_')) {
      const clusterIndex = parseInt(selectedSubstance.replace('cluster_', ''))
      if (clusterIndex >= 0 && clusterIndex < clusters.length) {
        return `Cluster ${clusterIndex + 1}`
      }
    }
    
    return 'Largest Cluster'
  }

  const currentClusterTitle = getCurrentClusterInfo()

  return (
    <div className="graph-preview">
      <h3 className="text-md font-semibold text-white mb-2">{currentClusterTitle}</h3>
      
      {previewData && (
        <div className="mb-2">
          {renderGraph()}
          
          {stats && (
            <div className="mt-2 text-xs text-gray-300">
              <div>Size: {stats.size} particles</div>
              <div>Bonds: {stats.bonds}</div>
              <div>Density: {stats.density}</div>
              <div>Avg Mass: {stats.avgMass}</div>
              <div>Total Energy: {stats.totalEnergy}</div>
              <div>Avg Energy: {stats.avgEnergy}</div>
            </div>
          )}
          
          {/* Export buttons */}
          <div className="export-buttons">
            <button
              onClick={() => exportAsImage('png')}
              className="export-btn png"
              title="Save as PNG"
            >
              <Download size={12} />
              PNG
            </button>
            <button
              onClick={() => exportAsImage('jpg')}
              className="export-btn jpg"
              title="Save as JPG"
            >
              <Download size={12} />
              JPG
            </button>
          </div>
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
