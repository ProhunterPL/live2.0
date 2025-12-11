import './SLACompliance.css'

interface SLAComplianceProps {
  title: string
  data: {
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
}

export default function SLACompliance({ title, data }: SLAComplianceProps) {
  const getStatusColor = (status: string) => {
    if (status === 'compliant') return 'operational'
    if (status === 'warning') return 'degraded'
    return 'down'
  }

  return (
    <div className="status-card">
      <h2>{title}</h2>
      <div className="status-indicator-wrapper">
        <span className={`status-indicator ${getStatusColor(data.status)}`}></span>
        <span>{data.status.toUpperCase()}</span>
      </div>
      
      <div className="sla-metrics">
        <div className="metric">
          <div className="metric-label">Uptime</div>
          <div className="metric-value">
            {data.metrics.uptime.toFixed(2)}% / {data.metrics.uptime_target}%
          </div>
        </div>
        
        <div className="metric">
          <div className="metric-label">Response Time (p95)</div>
          <div className="metric-value">
            {data.metrics.response_time_p95.toFixed(0)}ms / {data.metrics.response_time_p95_target}ms
          </div>
        </div>
        
        <div className="metric">
          <div className="metric-label">Error Rate</div>
          <div className="metric-value">
            {data.metrics.error_rate.toFixed(2)}% / {data.metrics.error_rate_target}%
          </div>
        </div>
      </div>
      
      {data.violations.length > 0 && (
        <div className="violations">
          <strong>Violations:</strong>
          <ul>
            {data.violations.map((v, i) => (
              <li key={i}>{v}</li>
            ))}
          </ul>
        </div>
      )}
      
      {data.warnings.length > 0 && (
        <div className="warnings">
          <strong>Warnings:</strong>
          <ul>
            {data.warnings.map((w, i) => (
              <li key={i}>{w}</li>
            ))}
          </ul>
        </div>
      )}
    </div>
  )
}
