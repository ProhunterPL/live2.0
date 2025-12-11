import { Metrics } from '../types/index'
import './StatusCard.css'

interface StatusCardProps {
  title: string
  metrics: Metrics
}

export default function StatusCard({ title, metrics }: StatusCardProps) {
  const getStatusColor = (uptime: number) => {
    if (uptime >= 99.9) return 'operational'
    if (uptime >= 95) return 'degraded'
    return 'down'
  }

  const statusColor = getStatusColor(metrics.uptime)

  return (
    <div className="status-card">
      <h2>{title}</h2>
      <div className="status-indicator-wrapper">
        <span className={`status-indicator ${statusColor}`}></span>
        <span>
          {metrics.uptime >= 99.9 ? 'Operational' : 
           metrics.uptime >= 95 ? 'Degraded' : 'Down'}
        </span>
      </div>
      
      <div className="metrics">
        <div className="metric">
          <div className="metric-label">Uptime</div>
          <div className="metric-value">{metrics.uptime.toFixed(2)}%</div>
        </div>
        
        <div className="metric">
          <div className="metric-label">Response Time (p95)</div>
          <div className="metric-value">
            {metrics.response_time.p95.toFixed(0)}ms
          </div>
        </div>
        
        <div className="metric">
          <div className="metric-label">Error Rate</div>
          <div className="metric-value">
            {metrics.error_rate.toFixed(2)}%
          </div>
        </div>
      </div>
      
      <div className="timestamp">
        Last updated: {new Date(metrics.timestamp).toLocaleString()}
      </div>
    </div>
  )
}
