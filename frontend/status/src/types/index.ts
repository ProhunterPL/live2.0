export interface Metrics {
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
