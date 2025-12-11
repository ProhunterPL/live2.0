import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts'
import { Metrics } from '../types/index'
import './MetricsChart.css'

interface MetricsChartProps {
  title: string
  metrics: Metrics
}

export default function MetricsChart({ title, metrics }: MetricsChartProps) {
  // For now, show current metrics as single point
  // In production, this would show historical data
  const data = [
    {
      time: new Date(metrics.timestamp).toLocaleTimeString(),
      p50: metrics.response_time.p50,
      p95: metrics.response_time.p95,
      p99: metrics.response_time.p99,
    }
  ]

  return (
    <div className="status-card">
      <h2>{title}</h2>
      <ResponsiveContainer width="100%" height={300}>
        <LineChart data={data}>
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis dataKey="time" />
          <YAxis label={{ value: 'ms', angle: -90, position: 'insideLeft' }} />
          <Tooltip />
          <Legend />
          <Line type="monotone" dataKey="p50" stroke="#8884d8" name="p50" />
          <Line type="monotone" dataKey="p95" stroke="#82ca9d" name="p95" />
          <Line type="monotone" dataKey="p99" stroke="#ffc658" name="p99" />
        </LineChart>
      </ResponsiveContainer>
    </div>
  )
}
