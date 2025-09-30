// Live 2.0 API Client

import type { 
  SimulationStatus, 
  SimulationData,
  CreateSimulationRequest,
  UpdateSimulationRequest 
} from './types'

export class SimulationAPI {
  private baseUrl: string

  constructor(baseUrl: string = 'http://localhost:8001') {
    this.baseUrl = baseUrl
  }

  async createSimulation(request: CreateSimulationRequest): Promise<{ success: boolean; simulation_id?: string; message?: string }> {
    const response = await fetch(`${this.baseUrl}/simulation/create`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        config: request.config,
        mode: request.config.mode
      }),
    })

    if (!response.ok) {
      throw new Error(`Failed to create simulation: ${response.statusText}`)
    }

    return response.json()
  }

  async getSimulationStatus(id: string): Promise<any> {
    const response = await fetch(`${this.baseUrl}/simulation/${id}/status`)

    if (!response.ok) {
      throw new Error(`Failed to get simulation: ${response.statusText}`)
    }
    return response.json()
  }

  async updateSimulation(id: string, request: UpdateSimulationRequest): Promise<any> {
    // Backend doesn't have a direct update endpoint, so we'll use the specific action endpoints
    if (request.status === 'running') {
      return this.startSimulation(id)
    } else if (request.status === 'paused') {
      return this.pauseSimulation(id)
    } else if (request.status === 'stopped') {
      return this.stopSimulation(id)
    }
    
    // If no status change, just return current status
    return this.getSimulationStatus(id)
  }

  async deleteSimulation(id: string): Promise<void> {
    // Backend doesn't have a delete endpoint, so we'll just stop the simulation
    await this.stopSimulation(id)
  }

  async startSimulation(id: string): Promise<any> {
    const response = await fetch(`${this.baseUrl}/simulation/${id}/start`, {
      method: 'POST',
    })

    if (!response.ok) {
      throw new Error(`Failed to start simulation: ${response.statusText}`)
    }

    return this.getSimulationStatus(id)
  }

  async pauseSimulation(id: string): Promise<any> {
    const response = await fetch(`${this.baseUrl}/simulation/${id}/pause`, {
      method: 'POST',
    })

    if (!response.ok) {
      throw new Error(`Failed to pause simulation: ${response.statusText}`)
    }

    return this.getSimulationStatus(id)
  }

  async stopSimulation(id: string): Promise<any> {
    const response = await fetch(`${this.baseUrl}/simulation/${id}/stop`, {
      method: 'POST',
    })

    if (!response.ok) {
      throw new Error(`Failed to stop simulation: ${response.statusText}`)
    }

    return this.getSimulationStatus(id)
  }

  async getSimulationData(_id: string): Promise<SimulationData> {
    // Backend doesn't have a direct data endpoint, so we'll return mock data for now
    // In a real implementation, this would come from WebSocket or a different endpoint
    return {
      // Provide empty structures matching expected shape
      particles: { positions: [], attributes: [], active_mask: [] },
      bonds: [],
      energy_field: [],
      metrics: {
        particle_count: 0,
        bond_count: 0,
        cluster_count: 0,
        novelty_rate: 0
      },
      timestamp: Date.now()
    }
  }

  async getSimulations(): Promise<SimulationStatus[]> {
    // Backend doesn't have a list endpoint, so we'll return empty array
    return []
  }

  async resumeSimulation(id: string): Promise<any> {
    // Backend has separate resume endpoint
    const response = await fetch(`${this.baseUrl}/simulation/${id}/resume`, { method: 'POST' })
    if (!response.ok) throw new Error(`Failed to resume simulation: ${response.statusText}`)
    return this.getSimulationStatus(id)
  }

  async resetSimulation(id: string): Promise<any> {
    const response = await fetch(`${this.baseUrl}/simulation/${id}/reset`, { method: 'POST' })
    if (!response.ok) throw new Error(`Failed to reset simulation: ${response.statusText}`)
    return this.getSimulationStatus(id)
  }

  async getHealth(): Promise<{ status: string; message: string }> {
    const response = await fetch(`${this.baseUrl}/`)

    if (!response.ok) {
      throw new Error(`Failed to get health: ${response.statusText}`)
    }

    return response.json()
  }

  async getMetrics(id: string): Promise<any> {
    const response = await fetch(`${this.baseUrl}/simulation/${id}/metrics`)
    if (!response.ok) throw new Error(`Failed to get metrics: ${response.statusText}`)
    return response.json()
  }

  async getNovelSubstances(id: string, count = 5): Promise<any> {
    const response = await fetch(`${this.baseUrl}/simulation/${id}/novel-substances?count=${count}`)
    if (!response.ok) throw new Error(`Failed to get novel substances: ${response.statusText}`)
    return response.json()
  }

  async saveSnapshot(id: string, filename: string): Promise<any> {
    const response = await fetch(`${this.baseUrl}/simulation/${id}/snapshot/save?filename=${encodeURIComponent(filename)}`, { method: 'POST' })
    if (!response.ok) throw new Error(`Failed to save snapshot: ${response.statusText}`)
    return response.json()
  }

  async loadSnapshot(id: string, filename: string): Promise<any> {
    const response = await fetch(`${this.baseUrl}/simulation/${id}/snapshot/load?filename=${encodeURIComponent(filename)}`, { method: 'POST' })
    if (!response.ok) throw new Error(`Failed to load snapshot: ${response.statusText}`)
    return response.json()
  }
}
