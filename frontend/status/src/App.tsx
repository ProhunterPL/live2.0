import { useState, useEffect } from 'react'
import axios from 'axios'
import StatusCard from './components/StatusCard'
import MetricsChart from './components/MetricsChart'
import SLACompliance from './components/SLACompliance'
import './App.css'

const API_BASE = import.meta.env.VITE_API_BASE || 'http://localhost:8000'

interface Metrics {
  project: string
  tier: string
  uptime: number
  response_time: {
    p50: number
    p95: number
    p99: number
  }
  error_rate: number
  timestamp: string
}

interface SLAComplianceData {
  status: string
  month: string
  project: string
  tier: string
  metrics: {
    uptime: number
    uptime_target: number
    response_time_p95: number
    response_time_p95_target: number
    error_rate: number
    error_rate_target: number
  }
  violations: string[]
  warnings: string[]
  compliant: boolean
}

function App() {
  const [live2Metrics, setLive2Metrics] = useState<Metrics | null>(null)
  const [legallyMetrics, setLegallyMetrics] = useState<Metrics | null>(null)
  const [live2SLA, setLive2SLA] = useState<SLAComplianceData | null>(null)
  const [legallySLA, setLegallySLA] = useState<SLAComplianceData | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  const fetchMetrics = async () => {
    try {
      setError(null)
      
      // Fetch Live 2.0 metrics
      const live2Response = await axios.get(`${API_BASE}/status/metrics`, {
        params: { project: 'live2', tier: 'hobby' }
      })
      setLive2Metrics(live2Response.data)
      
      // Fetch Legally metrics
      const legallyResponse = await axios.get(`${API_BASE}/status/metrics`, {
        params: { project: 'legally', tier: 'starter' }
      })
      setLegallyMetrics(legallyResponse.data)
      
      // Fetch SLA compliance
      const live2SLAResponse = await axios.get(`${API_BASE}/status/sla`, {
        params: { project: 'live2', tier: 'hobby' }
      })
      setLive2SLA(live2SLAResponse.data)
      
      const legallySLAResponse = await axios.get(`${API_BASE}/status/sla`, {
        params: { project: 'legally', tier: 'starter' }
      })
      setLegallySLA(legallySLAResponse.data)
      
      setLoading(false)
    } catch (err: any) {
      setError(err.message || 'Failed to fetch metrics')
      setLoading(false)
    }
  }

  useEffect(() => {
    fetchMetrics()
    // Refresh every 30 seconds
    const interval = setInterval(fetchMetrics, 30000)
    return () => clearInterval(interval)
  }, [])

  if (loading) {
    return (
      <div className="container">
        <div className="loading">Loading status...</div>
      </div>
    )
  }

  return (
    <div className="container">
      <div className="header">
        <h1>Live 2.0 & Legally - System Status</h1>
        <p>Real-time monitoring and SLA compliance</p>
      </div>

      {error && <div className="error">Error: {error}</div>}

      <div className="status-grid">
        {live2Metrics && (
          <StatusCard
            title="Live 2.0"
            metrics={live2Metrics}
          />
        )}
        
        {legallyMetrics && (
          <StatusCard
            title="Legally"
            metrics={legallyMetrics}
          />
        )}
      </div>

      <div className="status-grid">
        {live2SLA && (
          <SLACompliance
            title="Live 2.0 SLA"
            data={live2SLA}
          />
        )}
        
        {legallySLA && (
          <SLACompliance
            title="Legally SLA"
            data={legallySLA}
          />
        )}
      </div>

      {live2Metrics && (
        <MetricsChart
          title="Live 2.0 Response Times"
          metrics={live2Metrics}
        />
      )}
    </div>
  )
}

export default App
