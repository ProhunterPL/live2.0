/**
 * Test suite for Live 2.0 frontend components
 */

import React from 'react'
import { render, screen, fireEvent, waitFor } from '@testing-library/react'
import '@testing-library/jest-dom'
import App from '../App'
import HeatmapCanvas from '../components/HeatmapCanvas'
import Controls from '../components/Controls'
import GraphPreview from '../components/GraphPreview'
import { SimulationAPI } from '../lib/api'
import { WebSocketClient } from '../lib/ws'
import type { SimulationData, SimulationStatus } from '../lib/types'

// Mock the API and WebSocket clients
jest.mock('../lib/api')
jest.mock('../lib/ws')

const MockedSimulationAPI = SimulationAPI as jest.MockedClass<typeof SimulationAPI>
const MockedWebSocketClient = WebSocketClient as jest.MockedClass<typeof WebSocketClient>

describe('App Component', () => {
  beforeEach(() => {
    jest.clearAllMocks()
    
    // Mock API responses
    MockedSimulationAPI.prototype.createSimulation.mockResolvedValue({
      success: true,
      message: 'Simulation created successfully',
      simulation_id: 'sim_1234567890'
    })
    
    MockedSimulationAPI.prototype.getSimulationStatus.mockResolvedValue({
      simulation_id: 'sim_1234567890',
      is_running: true,
      is_paused: false,
      current_time: 45.2,
      step_count: 4520,
      particle_count: 150,
      novelty_rate: 0.15,
      health_score: 0.85
    })
    
    MockedSimulationAPI.prototype.startSimulation.mockResolvedValue({
      success: true,
      message: 'Simulation started'
    })
    
    MockedSimulationAPI.prototype.pauseSimulation.mockResolvedValue({
      success: true,
      message: 'Simulation paused'
    })
    
    MockedSimulationAPI.prototype.resumeSimulation.mockResolvedValue({
      success: true,
      message: 'Simulation resumed'
    })
    
    MockedSimulationAPI.prototype.stopSimulation.mockResolvedValue({
      success: true,
      message: 'Simulation stopped'
    })
    
    MockedSimulationAPI.prototype.resetSimulation.mockResolvedValue({
      success: true,
      message: 'Simulation reset'
    })
    
    // Mock WebSocket
    MockedWebSocketClient.prototype.connect.mockResolvedValue(undefined)
    MockedWebSocketClient.prototype.disconnect.mockImplementation(() => {})
    MockedWebSocketClient.prototype.isConnected = true
  })

  test('renders app without crashing', () => {
    render(<App />)
    expect(screen.getByText('Live 2.0 Simulation API')).toBeInTheDocument()
  })

  test('initializes simulation on mount', async () => {
    render(<App />)
    
    await waitFor(() => {
      expect(MockedSimulationAPI.prototype.createSimulation).toHaveBeenCalled()
    })
  })

  test('displays simulation status', async () => {
    render(<App />)
    
    await waitFor(() => {
      expect(screen.getByText('Running')).toBeInTheDocument()
    })
  })

  test('handles pause/resume button click', async () => {
    render(<App />)
    
    await waitFor(() => {
      const pauseButton = screen.getByText('Pause')
      fireEvent.click(pauseButton)
    })
    
    await waitFor(() => {
      expect(MockedSimulationAPI.prototype.pauseSimulation).toHaveBeenCalled()
    })
  })

  test('handles stop button click', async () => {
    render(<App />)
    
    await waitFor(() => {
      const stopButton = screen.getByText('Stop')
      fireEvent.click(stopButton)
    })
    
    await waitFor(() => {
      expect(MockedSimulationAPI.prototype.stopSimulation).toHaveBeenCalled()
    })
  })

  test('handles reset button click', async () => {
    render(<App />)
    
    await waitFor(() => {
      const resetButton = screen.getByText('Reset')
      fireEvent.click(resetButton)
    })
    
    await waitFor(() => {
      expect(MockedSimulationAPI.prototype.resetSimulation).toHaveBeenCalled()
    })
  })
})

describe('HeatmapCanvas Component', () => {
  const mockData: SimulationData = {
    particles: {
      positions: [[10, 20], [30, 40], [50, 60]],
      attributes: [[1.0, 0.5, -0.3, 0.1], [1.2, -0.2, 0.4, -0.1], [0.8, 0.1, 0.2, 0.0]],
      active_mask: [1, 1, 1]
    },
    energy_field: [
      [0.1, 0.2, 0.3],
      [0.4, 0.5, 0.6],
      [0.7, 0.8, 0.9]
    ],
    bonds: [
      { particle_i: 0, particle_j: 1, strength: 0.8 },
      { particle_i: 1, particle_j: 2, strength: 0.6 }
    ],
    clusters: [
      { particles: [0, 1, 2] }
    ],
    metrics: {
      particle_count: 3,
      total_energy: 100.0,
      total_mass: 3.0,
      bond_count: 2,
      cluster_count: 1,
      energy_field_sum: 4.5,
      energy_field_max: 0.9,
      energy_field_mean: 0.5,
      novelty_rate: 0.15,
      discovery_rate: 2.3,
      total_substances: 5,
      health_score: 0.85
    }
  }

  test('renders canvas element', () => {
    render(
      <HeatmapCanvas
        data={mockData}
        selectedSubstance={null}
        isConnected={true}
      />
    )
    
    const canvas = screen.getByRole('img', { hidden: true })
    expect(canvas).toBeInTheDocument()
  })

  test('displays connection status', () => {
    render(
      <HeatmapCanvas
        data={mockData}
        selectedSubstance={null}
        isConnected={true}
      />
    )
    
    expect(screen.getByText('Connected')).toBeInTheDocument()
  })

  test('displays disconnected status', () => {
    render(
      <HeatmapCanvas
        data={mockData}
        selectedSubstance={null}
        isConnected={false}
      />
    )
    
    expect(screen.getByText('Disconnected')).toBeInTheDocument()
  })

  test('renders control checkboxes', () => {
    render(
      <HeatmapCanvas
        data={mockData}
        selectedSubstance={null}
        isConnected={true}
      />
    )
    
    expect(screen.getByLabelText('Particles')).toBeInTheDocument()
    expect(screen.getByLabelText('Energy Field')).toBeInTheDocument()
    expect(screen.getByLabelText('Bonds')).toBeInTheDocument()
  })

  test('handles checkbox changes', () => {
    render(
      <HeatmapCanvas
        data={mockData}
        selectedSubstance={null}
        isConnected={true}
      />
    )
    
    const particlesCheckbox = screen.getByLabelText('Particles')
    fireEvent.click(particlesCheckbox)
    
    expect(particlesCheckbox).not.toBeChecked()
  })

  test('handles mouse movement', () => {
    render(
      <HeatmapCanvas
        data={mockData}
        selectedSubstance={null}
        isConnected={true}
      />
    )
    
    const canvas = screen.getByRole('img', { hidden: true })
    fireEvent.mouseMove(canvas, { clientX: 100, clientY: 100 })
    
    // Should not crash
    expect(canvas).toBeInTheDocument()
  })
})

describe('Controls Component', () => {
  const mockStatus: SimulationStatus = {
    simulation_id: 'sim_1234567890',
    is_running: true,
    is_paused: false,
    current_time: 45.2,
    step_count: 4520,
    particle_count: 150,
    novelty_rate: 0.15,
    health_score: 0.85
  }

  beforeEach(() => {
    MockedSimulationAPI.prototype.getMetrics.mockResolvedValue({
      metrics: {
        particle_count: 150,
        total_energy: 1250.5,
        total_mass: 180.3,
        bond_count: 75,
        cluster_count: 12,
        energy_field_sum: 5000.0,
        energy_field_max: 25.5,
        energy_field_mean: 0.076,
        novelty_rate: 0.15,
        discovery_rate: 2.3,
        total_substances: 45,
        health_score: 0.85
      }
    })
    
    MockedSimulationAPI.prototype.getNovelSubstances.mockResolvedValue({
      substances: [
        {
          id: 'SUB_abc12345_1234567890',
          timestamp: 45.2,
          size: 5,
          complexity: 12.3,
          properties: {
            density: 0.4,
            diameter: 3.2
          },
          graph: {
            particles: [1, 2, 3, 4, 5],
            bonds: [[1, 2], [2, 3], [3, 4], [4, 5]],
            particle_attributes: {
              1: [1.0, 0.5, -0.3, 0.1],
              2: [1.2, -0.2, 0.4, -0.1]
            }
          }
        }
      ]
    })
    
    MockedSimulationAPI.prototype.saveSnapshot.mockResolvedValue({
      success: true,
      filename: 'test_snapshot.json'
    })
    
    MockedSimulationAPI.prototype.loadSnapshot.mockResolvedValue({
      success: true,
      message: 'Snapshot loaded'
    })
  })

  test('renders controls without crashing', () => {
    render(
      <Controls
        simulationId="sim_1234567890"
        status={mockStatus}
        onStatusUpdate={jest.fn()}
      />
    )
    
    expect(screen.getByText('Simulation Controls')).toBeInTheDocument()
  })

  test('displays simulation status', () => {
    render(
      <Controls
        simulationId="sim_1234567890"
        status={mockStatus}
        onStatusUpdate={jest.fn()}
      />
    )
    
    expect(screen.getByText('Running')).toBeInTheDocument()
  })

  test('displays simulation time', () => {
    render(
      <Controls
        simulationId="sim_1234567890"
        status={mockStatus}
        onStatusUpdate={jest.fn()}
      />
    )
    
    expect(screen.getByText('0:45')).toBeInTheDocument()
  })

  test('displays step count', () => {
    render(
      <Controls
        simulationId="sim_1234567890"
        status={mockStatus}
        onStatusUpdate={jest.fn()}
      />
    )
    
    expect(screen.getByText('4.5K')).toBeInTheDocument()
  })

  test('displays metrics', async () => {
    render(
      <Controls
        simulationId="sim_1234567890"
        status={mockStatus}
        onStatusUpdate={jest.fn()}
      />
    )
    
    await waitFor(() => {
      expect(screen.getByText('150')).toBeInTheDocument() // particle count
      expect(screen.getByText('75')).toBeInTheDocument()  // bond count
      expect(screen.getByText('12')).toBeInTheDocument()  // cluster count
    })
  })

  test('handles snapshot save', async () => {
    render(
      <Controls
        simulationId="sim_1234567890"
        status={mockStatus}
        onStatusUpdate={jest.fn()}
      />
    )
    
    const filenameInput = screen.getByPlaceholderText('Snapshot filename')
    fireEvent.change(filenameInput, { target: { value: 'test_snapshot.json' } })
    
    const saveButton = screen.getByText('Save')
    fireEvent.click(saveButton)
    
    await waitFor(() => {
      expect(MockedSimulationAPI.prototype.saveSnapshot).toHaveBeenCalledWith(
        'sim_1234567890',
        'test_snapshot.json'
      )
    })
  })

  test('handles snapshot load', async () => {
    render(
      <Controls
        simulationId="sim_1234567890"
        status={mockStatus}
        onStatusUpdate={jest.fn()}
      />
    )
    
    const filenameInput = screen.getByPlaceholderText('Snapshot filename')
    fireEvent.change(filenameInput, { target: { value: 'test_snapshot.json' } })
    
    const loadButton = screen.getByText('Load')
    fireEvent.click(loadButton)
    
    await waitFor(() => {
      expect(MockedSimulationAPI.prototype.loadSnapshot).toHaveBeenCalledWith(
        'sim_1234567890',
        'test_snapshot.json'
      )
    })
  })

  test('toggles settings visibility', () => {
    render(
      <Controls
        simulationId="sim_1234567890"
        status={mockStatus}
        onStatusUpdate={jest.fn()}
      />
    )
    
    const settingsButton = screen.getByText('Show Settings')
    fireEvent.click(settingsButton)
    
    expect(screen.getByText('Hide Settings')).toBeInTheDocument()
    expect(screen.getByText('Advanced Settings')).toBeInTheDocument()
  })
})

describe('GraphPreview Component', () => {
  const mockData: SimulationData = {
    particles: {
      positions: [[10, 20], [30, 40], [50, 60], [70, 80]],
      attributes: [
        [1.0, 0.5, -0.3, 0.1],
        [1.2, -0.2, 0.4, -0.1],
        [0.8, 0.1, 0.2, 0.0],
        [1.1, -0.1, 0.3, 0.2]
      ],
      active_mask: [1, 1, 1, 1]
    },
    energy_field: [
      [0.1, 0.2, 0.3],
      [0.4, 0.5, 0.6],
      [0.7, 0.8, 0.9]
    ],
    bonds: [
      { particle_i: 0, particle_j: 1, strength: 0.8 },
      { particle_i: 1, particle_j: 2, strength: 0.6 },
      { particle_i: 2, particle_j: 3, strength: 0.7 },
      { particle_i: 3, particle_j: 0, strength: 0.5 }
    ],
    clusters: [
      { particles: [0, 1, 2, 3] }
    ],
    metrics: {
      particle_count: 4,
      total_energy: 100.0,
      total_mass: 4.1,
      bond_count: 4,
      cluster_count: 1,
      energy_field_sum: 4.5,
      energy_field_max: 0.9,
      energy_field_mean: 0.5,
      novelty_rate: 0.15,
      discovery_rate: 2.3,
      total_substances: 5,
      health_score: 0.85
    }
  }

  test('renders graph preview without crashing', () => {
    render(
      <GraphPreview
        data={mockData}
        selectedSubstance={null}
        onSubstanceSelect={jest.fn()}
      />
    )
    
    expect(screen.getByText('Largest Cluster')).toBeInTheDocument()
  })

  test('renders SVG graph', () => {
    render(
      <GraphPreview
        data={mockData}
        selectedSubstance={null}
        onSubstanceSelect={jest.fn()}
      />
    )
    
    const svg = screen.getByRole('img', { hidden: true })
    expect(svg).toBeInTheDocument()
  })

  test('displays cluster statistics', () => {
    render(
      <GraphPreview
        data={mockData}
        selectedSubstance={null}
        onSubstanceSelect={jest.fn()}
      />
    )
    
    expect(screen.getByText('Size: 4 particles')).toBeInTheDocument()
    expect(screen.getByText('Bonds: 4')).toBeInTheDocument()
  })

  test('handles cluster selection', () => {
    const onSubstanceSelect = jest.fn()
    
    render(
      <GraphPreview
        data={mockData}
        selectedSubstance={null}
        onSubstanceSelect={onSubstanceSelect}
      />
    )
    
    const clusterButton = screen.getByText('Cluster 1: 4 particles')
    fireEvent.click(clusterButton)
    
    expect(onSubstanceSelect).toHaveBeenCalledWith('cluster_0')
  })

  test('displays instructions', () => {
    render(
      <GraphPreview
        data={mockData}
        selectedSubstance={null}
        onSubstanceSelect={jest.fn()}
      />
    )
    
    expect(screen.getByText(/Hover over particles in the main view/)).toBeInTheDocument()
    expect(screen.getByText(/Click on clusters to highlight them/)).toBeInTheDocument()
  })

  test('handles empty data gracefully', () => {
    const emptyData: SimulationData = {
      particles: {
        positions: [],
        attributes: [],
        active_mask: []
      },
      energy_field: [],
      bonds: [],
      clusters: [],
      metrics: {
        particle_count: 0,
        total_energy: 0,
        total_mass: 0,
        bond_count: 0,
        cluster_count: 0,
        energy_field_sum: 0,
        energy_field_max: 0,
        energy_field_mean: 0,
        novelty_rate: 0,
        discovery_rate: 0,
        total_substances: 0,
        health_score: 0
      }
    }
    
    render(
      <GraphPreview
        data={emptyData}
        selectedSubstance={null}
        onSubstanceSelect={jest.fn()}
      />
    )
    
    expect(screen.getByText('No clusters to display')).toBeInTheDocument()
  })
})

describe('API Client', () => {
  test('creates simulation', async () => {
    const api = new SimulationAPI()
    
    MockedSimulationAPI.prototype.createSimulation.mockResolvedValue({
      success: true,
      message: 'Simulation created successfully',
      simulation_id: 'sim_1234567890'
    })
    
    const result = await api.createSimulation({
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
        seed: 42
      },
      mode: 'open_chemistry'
    })
    
    expect(result.success).toBe(true)
    expect(result.simulation_id).toBe('sim_1234567890')
  })

  test('handles API errors', async () => {
    const api = new SimulationAPI()
    
    MockedSimulationAPI.prototype.createSimulation.mockRejectedValue(
      new Error('API Error')
    )
    
    await expect(api.createSimulation({
      config: {} as any,
      mode: 'open_chemistry'
    })).rejects.toThrow('API Error')
  })
})

describe('WebSocket Client', () => {
  test('connects to WebSocket', async () => {
    const wsClient = new WebSocketClient()
    
    MockedWebSocketClient.prototype.connect.mockResolvedValue(undefined)
    
    await wsClient.connect('ws://localhost:8000/simulation/sim_123/stream')
    
    expect(MockedWebSocketClient.prototype.connect).toHaveBeenCalledWith(
      'ws://localhost:8000/simulation/sim_123/stream'
    )
  })

  test('disconnects from WebSocket', () => {
    const wsClient = new WebSocketClient()
    
    wsClient.disconnect()
    
    expect(MockedWebSocketClient.prototype.disconnect).toHaveBeenCalled()
  })

  test('handles connection errors', async () => {
    const wsClient = new WebSocketClient()
    
    MockedWebSocketClient.prototype.connect.mockRejectedValue(
      new Error('Connection failed')
    )
    
    await expect(wsClient.connect('ws://invalid-url')).rejects.toThrow(
      'Connection failed'
    )
  })
})
