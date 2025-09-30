import React, { useRef, useEffect, useState, useMemo } from 'react'
import type { SimulationData } from '../lib/types'

interface HeatmapCanvasProps {
  data: SimulationData | null
  selectedSubstance: string | null
  isConnected: boolean
}

const HeatmapCanvas: React.FC<HeatmapCanvasProps> = ({
  data,
  selectedSubstance,
  isConnected
}) => {
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const [mousePos, setMousePos] = useState<{ x: number; y: number } | null>(null)
  const [showParticles, setShowParticles] = useState(true)
  const [showEnergy, setShowEnergy] = useState(true)
  const [showBonds, setShowBonds] = useState(true)
  const [particleSize, setParticleSize] = useState(2)
  const [energyOpacity, setEnergyOpacity] = useState(0.7)

  // Memoize expensive calculations
  const memoizedData = useMemo(() => {
    if (!data) return null
    
    // Downsampling disabled for debugging

    return {
      ...data,
      energy_field: data.energy_field, // Disabled downsampling for debugging
      concentrations: data.concentrations
    }
  }, [data])

  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return

    const ctx = canvas.getContext('2d')
    if (!ctx) return

    // Set canvas size
    const resizeCanvas = () => {
      const rect = canvas.getBoundingClientRect()
      // Use 1x pixel ratio in preview/prod to reduce memory and ensure UI overlays visible
      const ratio = 1
      canvas.width = rect.width * ratio
      canvas.height = rect.height * ratio
      ctx.setTransform(1, 0, 0, 1, 0, 0)
      ctx.scale(ratio, ratio)
    }

    resizeCanvas()
    window.addEventListener('resize', resizeCanvas)

    return () => {
      window.removeEventListener('resize', resizeCanvas)
    }
  }, [])

  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas || !memoizedData) {
      console.log('HeatmapCanvas: No canvas or data', { canvas: !!canvas, data: !!memoizedData })
      return
    }

    console.log('HeatmapCanvas: Rendering with data:', {
      particles: memoizedData.particles,
      energy_field: memoizedData.energy_field?.length,
      bonds: memoizedData.bonds?.length,
      showParticles,
      showEnergy,
      showBonds
    })

    const ctx = canvas.getContext('2d')
    if (!ctx) return

    const rect = canvas.getBoundingClientRect()
    const width = rect.width
    const height = rect.height

    // Clear canvas
    ctx.fillStyle = '#000000'
    ctx.fillRect(0, 0, width, height)

    // Draw energy field
    if (showEnergy && memoizedData.energy_field) {
      drawEnergyField(ctx, memoizedData.energy_field, width, height)
    }

    // Draw concentration field for preset mode
    if (memoizedData.concentrations) {
      drawConcentrationField(ctx, memoizedData.concentrations, selectedSubstance, width, height)
    }

    // Draw bonds
    if (showBonds && memoizedData.bonds && memoizedData.particles) {
      let bonds: Array<{ particle_i: number; particle_j: number; strength: number }>
      
      if (Array.isArray(memoizedData.bonds) && memoizedData.bonds.length > 0) {
        if (Array.isArray(memoizedData.bonds[0])) {
          // Handle [number, number, number][] format
          bonds = (memoizedData.bonds as [number, number, number][]).map(([i, j, s]) => ({ 
            particle_i: i, 
            particle_j: j, 
            strength: s 
          }))
        } else {
          // Handle object format
          bonds = memoizedData.bonds as Array<{ particle_i: number; particle_j: number; strength: number }>
        }
      } else {
        bonds = []
      }
      
      drawBonds(ctx, bonds, memoizedData.particles.positions, width, height)
    }

    // Draw particles
    if (showParticles && memoizedData.particles) {
      drawParticles(ctx, memoizedData.particles, width, height)
    }

    // Draw mouse position info
    if (mousePos) {
      drawMouseInfo(ctx, mousePos, memoizedData, width, height)
    }

  }, [memoizedData, selectedSubstance, showParticles, showEnergy, showBonds, particleSize, energyOpacity, mousePos])

  // Throttle mouse position updates to reduce re-renders
  const throttledMouseMove = useRef<number | null>(null)
  const handleMouseMoveThrottled = (event: React.MouseEvent<HTMLCanvasElement>) => {
    if (throttledMouseMove.current) {
      clearTimeout(throttledMouseMove.current)
    }
    throttledMouseMove.current = window.setTimeout(() => {
      handleMouseMove(event)
    }, 16) // ~60 FPS
  }

  const drawEnergyField = (
    ctx: CanvasRenderingContext2D,
    energyField: number[][],
    width: number,
    height: number
  ) => {
    const fieldWidth = energyField[0]?.length || 0
    const fieldHeight = energyField.length || 0
    
    console.log('drawEnergyField:', { fieldWidth, fieldHeight, energyField: energyField.slice(0, 2) })
    console.log('Raw energy data sample:', energyField[0]?.slice(0, 10))
    
    if (fieldWidth === 0 || fieldHeight === 0) {
      console.log('drawEnergyField: Empty field')
      return
    }

    const cellWidth = width / fieldWidth
    const cellHeight = height / fieldHeight

    // Find max energy for normalization
    let maxEnergy = 0
    for (let y = 0; y < fieldHeight; y++) {
      for (let x = 0; x < fieldWidth; x++) {
        maxEnergy = Math.max(maxEnergy, energyField[y][x])
      }
    }

    console.log('drawEnergyField: maxEnergy =', maxEnergy)

    if (maxEnergy === 0) {
      console.log('drawEnergyField: No energy to draw')
      return
    }

    // Draw energy field
    for (let y = 0; y < fieldHeight; y++) {
      for (let x = 0; x < fieldWidth; x++) {
        const energy = energyField[y][x]
        const intensity = energy / maxEnergy
        
        if (intensity > 0) {
          // Color from blue (low) to red (high)
          const hue = (1 - intensity) * 240 // 240 = blue, 0 = red
          ctx.fillStyle = `hsla(${hue}, 100%, 50%, ${energyOpacity})`
          ctx.fillRect(x * cellWidth, y * cellHeight, cellWidth, cellHeight)
        }
      }
    }
  }

  const drawConcentrationField = (
    ctx: CanvasRenderingContext2D,
    concentrations: Record<string, number[][]>,
    species: string | null,
    width: number,
    height: number
  ) => {
    const keys = Object.keys(concentrations)
    if (keys.length === 0) return
    const key = species && concentrations[species] ? species : keys[0]
    const field = concentrations[key]
    const fieldHeight = field.length || 0
    const fieldWidth = field[0]?.length || 0
    if (fieldWidth === 0 || fieldHeight === 0) return

    const cellWidth = width / fieldWidth
    const cellHeight = height / fieldHeight

    // Find max concentration for normalization
    let maxVal = 0
    for (let y = 0; y < fieldHeight; y++) {
      for (let x = 0; x < fieldWidth; x++) {
        maxVal = Math.max(maxVal, field[y][x])
      }
    }
    if (maxVal === 0) return

    // Draw concentration heatmap (green palette)
    for (let y = 0; y < fieldHeight; y++) {
      for (let x = 0; x < fieldWidth; x++) {
        const v = field[y][x]
        const intensity = v / maxVal
        if (intensity > 0) {
          const hue = 120
          const light = 20 + 60 * intensity
          ctx.fillStyle = `hsla(${hue}, 100%, ${light}%, ${energyOpacity})`
          ctx.fillRect(x * cellWidth, y * cellHeight, cellWidth, cellHeight)
        }
      }
    }
  }

  const drawBonds = (
    ctx: CanvasRenderingContext2D,
    bonds: Array<{ particle_i: number; particle_j: number; strength: number }>,
    positions: [number, number][],
    width: number,
    height: number
  ) => {
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.3)'
    ctx.lineWidth = 1

    bonds.forEach(bond => {
      const pos1 = positions[bond.particle_i]
      const pos2 = positions[bond.particle_j]
      
      if (pos1 && pos2) {
        // Scale positions to canvas
        const x1 = (pos1[0] / 256) * width
        const y1 = (pos1[1] / 256) * height
        const x2 = (pos2[0] / 256) * width
        const y2 = (pos2[1] / 256) * height

        ctx.beginPath()
        ctx.moveTo(x1, y1)
        ctx.lineTo(x2, y2)
        ctx.stroke()
      }
    })
  }

  const drawParticles = (
    ctx: CanvasRenderingContext2D,
    particles: { positions: [number, number][]; attributes: [number, number, number, number][] },
    width: number,
    height: number
  ) => {
    console.log('drawParticles:', { 
      positions: particles.positions?.length, 
      attributes: particles.attributes?.length,
      samplePos: particles.positions?.slice(0, 3)
    })
    
    if (!particles.positions || !particles.attributes) {
      console.log('drawParticles: Missing positions or attributes')
      return
    }
    
    particles.positions.forEach((pos, index) => {
      const attr = particles.attributes[index]
      if (!attr || !pos) return

      const [mass, chargeX, chargeY, chargeZ] = attr
      
      // Scale position to canvas
      const x = (pos[0] / 256) * width
      const y = (pos[1] / 256) * height

      // Color based on charge
      const chargeMagnitude = Math.sqrt(chargeX * chargeX + chargeY * chargeY + chargeZ * chargeZ)
      const hue = chargeMagnitude > 0 ? (chargeX > 0 ? 0 : 120) : 240 // Red for positive, green for negative, blue for neutral
      
      ctx.fillStyle = `hsl(${hue}, 70%, 60%)`
      ctx.beginPath()
      ctx.arc(x, y, particleSize, 0, 2 * Math.PI)
      ctx.fill()

      // Draw mass indicator (size)
      const massSize = Math.max(1, mass * 2)
      ctx.strokeStyle = `hsl(${hue}, 70%, 80%)`
      ctx.lineWidth = 0.5
      ctx.beginPath()
      ctx.arc(x, y, massSize, 0, 2 * Math.PI)
      ctx.stroke()
    })
  }

  const drawMouseInfo = (
    ctx: CanvasRenderingContext2D,
    mousePos: { x: number; y: number },
    data: SimulationData,
    width: number,
    height: number
  ) => {
    // Convert mouse position to simulation coordinates
    const simX = (mousePos.x / width) * 256
    const simY = (mousePos.y / height) * 256

    // Find nearest particle
    let nearestParticle: { index: number; pos: [number, number]; attr: [number, number, number, number] } | null = null
    let minDistance = Infinity

    const parts = data.particles
    if (parts && parts.positions && parts.attributes) {
      parts.positions.forEach((pos, index) => {
        if (!pos) return
        const distance = Math.sqrt((pos[0] - simX) ** 2 + (pos[1] - simY) ** 2)
        if (distance < minDistance) {
          minDistance = distance
          nearestParticle = { index, pos, attr: parts.attributes[index] }
        }
      })
    }

    // Draw info box
    if (nearestParticle && minDistance < 10) {
      const p = nearestParticle as { index: number; pos: [number, number]; attr: [number, number, number, number] }
      const [mass, chargeX, chargeY, chargeZ] = p.attr

      ctx.fillStyle = 'rgba(0, 0, 0, 0.8)'
      ctx.fillRect(mousePos.x + 10, mousePos.y - 60, 150, 50)

      ctx.fillStyle = '#ffffff'
      ctx.font = '12px monospace'
      ctx.fillText(`Mass: ${mass.toFixed(2)}`, mousePos.x + 15, mousePos.y - 40)
      ctx.fillText(`Charge: (${chargeX.toFixed(2)}, ${chargeY.toFixed(2)}, ${chargeZ.toFixed(2)})`, mousePos.x + 15, mousePos.y - 25)
      ctx.fillText(`Pos: (${p.pos[0].toFixed(1)}, ${p.pos[1].toFixed(1)})`, mousePos.x + 15, mousePos.y - 10)
    }
  }

  const handleMouseMove = (event: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current
    if (!canvas) return

    const rect = canvas.getBoundingClientRect()
    setMousePos({
      x: event.clientX - rect.left,
      y: event.clientY - rect.top
    })
  }

  const handleMouseLeave = () => {
    setMousePos(null)
  }

  return (
    <div className="relative w-full h-full">
      <canvas
        ref={canvasRef}
        className="heatmap-canvas"
        onMouseMove={handleMouseMoveThrottled}
        onMouseLeave={handleMouseLeave}
      />
      
      {/* Controls overlay */}
      <div className="absolute bg-black bg-opacity-50 p-2 rounded" style={{ zIndex: 1000, top: 8, left: 8, position: 'absolute' }}>
        <div className="flex flex-col gap-2 text-sm">
          <label className="flex items-center gap-2">
            <input
              type="checkbox"
              checked={showParticles}
              onChange={(e) => setShowParticles(e.target.checked)}
            />
            Particles
          </label>
          <label className="flex items-center gap-2">
            <input
              type="checkbox"
              checked={showEnergy}
              onChange={(e) => setShowEnergy(e.target.checked)}
            />
            Energy Field
          </label>
          <label className="flex items-center gap-2">
            <input
              type="checkbox"
              checked={showBonds}
              onChange={(e) => setShowBonds(e.target.checked)}
            />
            Bonds
          </label>
          <div className="flex items-center gap-2">
            <label>Particle Size:</label>
            <input
              type="range"
              min="1"
              max="5"
              value={particleSize}
              onChange={(e) => setParticleSize(Number(e.target.value))}
              className="w-16"
            />
          </div>
          <div className="flex items-center gap-2">
            <label>Energy Opacity:</label>
            <input
              type="range"
              min="0"
              max="1"
              step="0.1"
              value={energyOpacity}
              onChange={(e) => setEnergyOpacity(Number(e.target.value))}
              className="w-16"
            />
          </div>
        </div>
      </div>

      {/* Connection status */}
      <div className="absolute top-4 right-4">
        <div className={`px-3 py-1 rounded text-sm ${
          isConnected ? 'bg-green-500 bg-opacity-20 text-green-400' : 'bg-red-500 bg-opacity-20 text-red-400'
        }`}>
          {isConnected ? 'Connected' : 'Disconnected'}
        </div>
      </div>
    </div>
  )
}

export default HeatmapCanvas
