/**
 * Component for displaying and managing API v1 jobs
 * Shows job status, progress, and download links for results stored in Supabase Storage
 */

import React, { useState, useEffect } from 'react'
import { APIv1Client, JobStatus } from '../lib/api_v1'
import { Download, RefreshCw, Clock, CheckCircle, XCircle, Loader, LogIn, UserPlus } from 'lucide-react'

interface APIv1JobsProps {
  apiBaseUrl?: string
  apiKey?: string
}

interface AuthResponse {
  access_token: string
  api_key: string
  user: {
    id: string
    email: string
    tier: string
    subscription_status: string
    api_key: string
  }
}

const APIv1Jobs: React.FC<APIv1JobsProps> = ({ 
  apiBaseUrl = 'http://localhost:8001',
  apiKey 
}) => {
  const [client, setClient] = useState<APIv1Client | null>(null)
  const [jobs, setJobs] = useState<JobStatus[]>([])
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [apiKeyInput, setApiKeyInput] = useState(apiKey || '')
  const [showAuth, setShowAuth] = useState(false)
  const [authMode, setAuthMode] = useState<'login' | 'register'>('login')
  const [email, setEmail] = useState('')
  const [password, setPassword] = useState('')
  const [tier, setTier] = useState('hobby')
  const [authLoading, setAuthLoading] = useState(false)
  const [user, setUser] = useState<AuthResponse['user'] | null>(null)

  useEffect(() => {
    const apiClient = new APIv1Client(apiBaseUrl, apiKey)
    setClient(apiClient)
    if (apiKey) {
      apiClient.setApiKey(apiKey)
    }
  }, [apiBaseUrl, apiKey])

  useEffect(() => {
    if (client) {
      loadJobs()
      // Poll for updates every 5 seconds
      const interval = setInterval(() => {
        loadJobs()
      }, 5000)
      
      return () => {
        clearInterval(interval)
      }
    }
  }, [client])

  const loadJobs = async () => {
    if (!client) return
    
    setLoading(true)
    setError(null)
    
    try {
      const jobList = await client.listJobs(20)
      setJobs(jobList)
    } catch (err: any) {
      // If endpoint doesn't exist, show empty list
      if (err.message?.includes('404') || err.message?.includes('not found')) {
        setJobs([])
      } else {
        setError(err.message || 'Failed to load jobs')
      }
    } finally {
      setLoading(false)
    }
  }

  const handleApiKeySubmit = (e: React.FormEvent) => {
    e.preventDefault()
    if (client && apiKeyInput.trim()) {
      client.setApiKey(apiKeyInput.trim())
      loadJobs()
    }
  }

  const handleRegister = async (e: React.FormEvent) => {
    e.preventDefault()
    setAuthLoading(true)
    setError(null)
    
    try {
      const response = await fetch(`${apiBaseUrl}/api/v1/auth/register`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ email, password, tier })
      })
      
      if (!response.ok) {
        const error = await response.json().catch(() => ({ detail: response.statusText }))
        throw new Error(error.detail || 'Registration failed')
      }
      
      const data: AuthResponse = await response.json()
      setUser(data.user)
      setApiKeyInput(data.api_key)
      if (client) {
        client.setApiKey(data.api_key)
      }
      setShowAuth(false)
      setEmail('')
      setPassword('')
      loadJobs()
    } catch (err: any) {
      setError(err.message || 'Registration failed')
    } finally {
      setAuthLoading(false)
    }
  }

  const handleLogin = async (e: React.FormEvent) => {
    e.preventDefault()
    setAuthLoading(true)
    setError(null)
    
    try {
      const response = await fetch(`${apiBaseUrl}/api/v1/auth/login`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ email, password })
      })
      
      if (!response.ok) {
        const error = await response.json().catch(() => ({ detail: response.statusText }))
        throw new Error(error.detail || 'Login failed')
      }
      
      const data: AuthResponse = await response.json()
      setUser(data.user)
      setApiKeyInput(data.api_key)
      if (client) {
        client.setApiKey(data.api_key)
      }
      setShowAuth(false)
      setEmail('')
      setPassword('')
      loadJobs()
    } catch (err: any) {
      setError(err.message || 'Login failed')
    } finally {
      setAuthLoading(false)
    }
  }

  const handleDownload = async (job: JobStatus) => {
    if (!client || !job.result_url) return
    
    try {
      const filename = `job_${job.job_id}_result.parquet`
      await client.downloadResult(job.result_url, filename)
    } catch (err: any) {
      alert(`Failed to download: ${err.message}`)
    }
  }

  const formatDate = (dateString: string) => {
    if (!dateString) return 'N/A'
    try {
      return new Date(dateString).toLocaleString()
    } catch {
      return dateString
    }
  }

  const getStatusIcon = (status: JobStatus['status']) => {
    switch (status) {
      case 'completed':
        return <CheckCircle className="w-5 h-5 text-green-500" />
      case 'failed':
        return <XCircle className="w-5 h-5 text-red-500" />
      case 'running':
        return <Loader className="w-5 h-5 text-blue-500 animate-spin" />
      default:
        return <Clock className="w-5 h-5 text-gray-500" />
    }
  }

  const getStatusColor = (status: JobStatus['status']) => {
    switch (status) {
      case 'completed':
        return 'bg-green-100 text-green-800'
      case 'failed':
        return 'bg-red-100 text-red-800'
      case 'running':
        return 'bg-blue-100 text-blue-800'
      default:
        return 'bg-gray-100 text-gray-800'
    }
  }

  if (!client) {
    return <div className="p-4">Initializing API client...</div>
  }

  return (
    <div className="p-6 max-w-6xl mx-auto">
      <div className="mb-6">
        <h2 className="text-2xl font-bold mb-4">API v1 Jobs</h2>
        <p className="text-gray-600 mb-4">
          View and manage your simulation jobs. Results are stored in Supabase Storage (S3-compatible).
        </p>
        
        {/* User Info */}
        {user && (
          <div className="mb-4 p-4 bg-green-50 border border-green-200 rounded-md">
            <p className="text-sm text-gray-700">
              <strong>Logged in as:</strong> {user.email} | <strong>Tier:</strong> {user.tier} | <strong>Status:</strong> {user.subscription_status}
            </p>
            <p className="text-xs text-gray-600 mt-1">
              <strong>API Key:</strong> {user.api_key}
            </p>
          </div>
        )}
        
        {/* Auth Buttons */}
        {!user && (
          <div className="mb-4 flex gap-2">
            <button
              onClick={() => {
                setShowAuth(true)
                setAuthMode('register')
              }}
              className="flex items-center gap-2 px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700"
            >
              <UserPlus className="w-4 h-4" />
              Register
            </button>
            <button
              onClick={() => {
                setShowAuth(true)
                setAuthMode('login')
              }}
              className="flex items-center gap-2 px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
            >
              <LogIn className="w-4 h-4" />
              Login
            </button>
          </div>
        )}
        
        {/* Auth Form */}
        {showAuth && (
          <div className="mb-4 p-4 border border-gray-300 rounded-md bg-gray-50">
            <h3 className="text-lg font-semibold mb-3">
              {authMode === 'register' ? 'Register New Account' : 'Login'}
            </h3>
            <form onSubmit={authMode === 'register' ? handleRegister : handleLogin}>
              <div className="space-y-3">
                <input
                  type="email"
                  value={email}
                  onChange={(e) => setEmail(e.target.value)}
                  placeholder="Email"
                  required
                  className="w-full px-4 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                />
                <input
                  type="password"
                  value={password}
                  onChange={(e) => setPassword(e.target.value)}
                  placeholder="Password"
                  required
                  className="w-full px-4 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                />
                {authMode === 'register' && (
                  <select
                    value={tier}
                    onChange={(e) => setTier(e.target.value)}
                    className="w-full px-4 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                  >
                    <option value="hobby">Hobby</option>
                    <option value="research">Research</option>
                    <option value="pro">Pro</option>
                  </select>
                )}
                <div className="flex gap-2">
                  <button
                    type="submit"
                    disabled={authLoading}
                    className="flex-1 px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 disabled:opacity-50"
                  >
                    {authLoading ? 'Loading...' : (authMode === 'register' ? 'Register' : 'Login')}
                  </button>
                  <button
                    type="button"
                    onClick={() => {
                      setShowAuth(false)
                      setError(null)
                    }}
                    className="px-4 py-2 bg-gray-300 text-gray-700 rounded-md hover:bg-gray-400"
                  >
                    Cancel
                  </button>
                </div>
              </div>
            </form>
          </div>
        )}
        
        {/* API Key Input (Manual) */}
        {!user && (
          <form onSubmit={handleApiKeySubmit} className="mb-4">
            <p className="text-sm text-gray-600 mb-2">Or enter API key manually:</p>
            <div className="flex gap-2">
              <input
                type="password"
                value={apiKeyInput}
                onChange={(e) => setApiKeyInput(e.target.value)}
                placeholder="Enter API Key (sk_live_...)"
                className="flex-1 px-4 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
              />
              <button
                type="submit"
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-500"
              >
                Set API Key
              </button>
            </div>
          </form>
        )}

        {/* Refresh Button */}
        <button
          onClick={loadJobs}
          disabled={loading}
          className="flex items-center gap-2 px-4 py-2 bg-gray-600 text-white rounded-md hover:bg-gray-700 disabled:opacity-50"
        >
          <RefreshCw className={`w-4 h-4 ${loading ? 'animate-spin' : ''}`} />
          Refresh
        </button>
      </div>

      {/* Error Message */}
      {error && (
        <div className="mb-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded">
          {error}
        </div>
      )}

      {/* Jobs List */}
      {jobs.length === 0 ? (
        <div className="p-8 text-center text-gray-500">
          {loading ? 'Loading jobs...' : 'No jobs found. Submit a job via API to see it here.'}
        </div>
      ) : (
        <div className="space-y-4">
          {jobs.map((job) => (
            <div
              key={job.job_id}
              className="border border-gray-300 rounded-lg p-4 hover:shadow-md transition-shadow"
            >
              <div className="flex items-start justify-between mb-3">
                <div className="flex items-center gap-3">
                  {getStatusIcon(job.status)}
                  <div>
                    <h3 className="font-semibold text-lg">Job {job.job_id.slice(0, 8)}...</h3>
                    <span className={`inline-block px-2 py-1 rounded text-sm font-medium ${getStatusColor(job.status)}`}>
                      {job.status.toUpperCase()}
                    </span>
                  </div>
                </div>
                
                {job.status === 'completed' && job.result_url && (
                  <button
                    onClick={() => handleDownload(job)}
                    className="flex items-center gap-2 px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700"
                  >
                    <Download className="w-4 h-4" />
                    Download Result
                  </button>
                )}
              </div>

              {/* Progress Bar */}
              {job.status === 'running' && (
                <div className="mb-3">
                  <div className="flex justify-between text-sm text-gray-600 mb-1">
                    <span>Progress</span>
                    <span>{job.progress}%</span>
                  </div>
                  <div className="w-full bg-gray-200 rounded-full h-2">
                    <div
                      className="bg-blue-600 h-2 rounded-full transition-all duration-300"
                      style={{ width: `${job.progress}%` }}
                    />
                  </div>
                </div>
              )}

              {/* Error Message */}
              {job.error && (
                <div className="mb-3 p-3 bg-red-50 border border-red-200 rounded text-sm text-red-700">
                  <strong>Error:</strong> {job.error}
                </div>
              )}

              {/* Result URL Info */}
              {job.result_url && (
                <div className="mb-3 p-3 bg-green-50 border border-green-200 rounded text-sm">
                  <strong>Result URL:</strong>{' '}
                  <a
                    href={job.result_url}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 hover:underline break-all"
                  >
                    {job.result_url.length > 80 ? `${job.result_url.slice(0, 80)}...` : job.result_url}
                  </a>
                  <p className="text-gray-600 mt-1 text-xs">
                    Stored in Supabase Storage (S3-compatible bucket: live-dysk)
                  </p>
                </div>
              )}

              {/* Metadata */}
              <div className="flex gap-4 text-sm text-gray-600">
                <div>
                  <strong>Created:</strong> {formatDate(job.created_at)}
                </div>
                {job.completed_at && (
                  <div>
                    <strong>Completed:</strong> {formatDate(job.completed_at)}
                  </div>
                )}
              </div>
            </div>
          ))}
        </div>
      )}
    </div>
  )
}

export default APIv1Jobs

