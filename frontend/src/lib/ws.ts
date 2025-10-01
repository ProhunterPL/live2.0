// Live 2.0 WebSocket Client

import type { SimulationData } from './types'
// eslint-disable-next-line @typescript-eslint/ban-ts-comment
// @ts-ignore - no types for msgpack-lite
import msgpack from 'msgpack-lite'

export class WebSocketClient {
  private ws: WebSocket | null = null
  private url: string
  private reconnectAttempts = 0
  private maxReconnectAttempts = 3  // Reduced from 5 to 3
  private reconnectDelay = 1000
  private listeners: Map<string, Function[]> = new Map()
  private shouldReconnect = true  // Flag to prevent reconnection for non-existent simulations

  constructor(url: string = 'ws://localhost:8001/simulation') {
    this.url = url
  }

  connect(urlOrId?: string): Promise<void> {
    return new Promise((resolve, reject) => {
      try {
        // Allow passing a full ws URL or just a simulationId
        const isFull = !!urlOrId && /^(ws|wss):\/\//.test(urlOrId)
        const wsUrl = isFull
          ? (urlOrId as string)
          : (urlOrId ? `${this.url}/${urlOrId}/stream` : this.url)
        
        // Reset reconnect flag for new connection attempts
        this.shouldReconnect = true
        
        console.log('WS connecting to:', wsUrl)
        this.ws = new WebSocket(wsUrl)
        this.ws.binaryType = 'arraybuffer'

        this.ws.onopen = () => {
          console.log('WebSocket connected')
          this.reconnectAttempts = 0
          this.emit('connected')
          resolve()
        }

        this.ws.onmessage = async (event) => {
          try {
            let buffer: ArrayBuffer
            if (event.data instanceof ArrayBuffer) {
              buffer = event.data
            } else if (event.data instanceof Blob) {
              buffer = await event.data.arrayBuffer()
            } else if (typeof event.data === 'string') {
              // Fallback: try JSON if server ever sends textual messages
              const json = JSON.parse(event.data)
              this.emit('simulation_data', json as SimulationData)
              return
            } else {
              console.warn('Unknown WebSocket message type')
              return
            }

            const data = msgpack.decode(new Uint8Array(buffer)) as SimulationData
            this.emit('simulation_data', data)
          } catch (error) {
            console.error('Failed to decode WebSocket message:', error)
          }
        }

        this.ws.onclose = (event) => {
          console.log('WebSocket disconnected:', event.code, event.reason)
          this.emit('disconnected')
          
          // Codes that indicate simulation doesn't exist or permanent failure:
          // 1008 = Policy Violation (backend returns this for non-existent simulation)
          // 1003 = Unsupported Data
          // 1002 = Protocol Error
          const permanentFailureCodes = [1008, 1003, 1002]
          
          if (permanentFailureCodes.includes(event.code)) {
            console.warn(`WebSocket closed with code ${event.code} (simulation not found or permanent error). Stopping reconnection attempts.`)
            this.shouldReconnect = false
            this.emit('simulation_not_found', { code: event.code, reason: event.reason })
            return
          }
          
          // Only reconnect if:
          // - it's not a normal closure (1000)
          // - we haven't exceeded max attempts
          // - we have a valid simulation ID
          // - shouldReconnect flag is true
          if (event.code !== 1000 && 
              this.shouldReconnect && 
              this.reconnectAttempts < this.maxReconnectAttempts && 
              urlOrId) {
            this.reconnect(urlOrId)
          } else if (this.reconnectAttempts >= this.maxReconnectAttempts) {
            console.warn(`Max reconnection attempts (${this.maxReconnectAttempts}) reached. Giving up.`)
            this.emit('reconnect_failed')
          }
        }

        this.ws.onerror = (error) => {
          console.error('WebSocket error:', error)
          this.emit('error', error)
          reject(error)
        }
      } catch (error) {
        reject(error)
      }
    })
  }

  private reconnect(simulationId?: string) {
    this.reconnectAttempts++
    console.log(`Attempting to reconnect (${this.reconnectAttempts}/${this.maxReconnectAttempts})...`)
    
    // Don't reconnect if simulationId is not provided
    if (!simulationId) {
      console.log('No simulation ID provided, skipping reconnect')
      return
    }
    
    setTimeout(() => {
      this.connect(simulationId).catch(console.error)
    }, this.reconnectDelay * this.reconnectAttempts)
  }

  // no-op legacy handler removed; messages handled in onmessage

  disconnect() {
    // Prevent automatic reconnection on manual disconnect
    this.shouldReconnect = false
    this.reconnectAttempts = 0
    
    if (this.ws) {
      this.ws.close(1000, 'Client disconnect')
      this.ws = null
    }
  }
  
  // Method to stop reconnection attempts without closing connection
  stopReconnecting() {
    this.shouldReconnect = false
    console.log('Reconnection attempts stopped')
  }

  send(message: any) {
    if (this.ws && this.ws.readyState === WebSocket.OPEN) {
      this.ws.send(JSON.stringify(message))
    } else {
      console.warn('WebSocket is not connected')
    }
  }

  on(event: string, callback: Function) {
    if (!this.listeners.has(event)) {
      this.listeners.set(event, [])
    }
    this.listeners.get(event)!.push(callback)
  }

  off(event: string, callback: Function) {
    const callbacks = this.listeners.get(event)
    if (callbacks) {
      const index = callbacks.indexOf(callback)
      if (index > -1) {
        callbacks.splice(index, 1)
      }
    }
  }

  private emit(event: string, data?: any) {
    const callbacks = this.listeners.get(event)
    if (callbacks) {
      callbacks.forEach(callback => callback(data))
    }
  }

  get isConnected(): boolean {
    return this.ws !== null && this.ws.readyState === WebSocket.OPEN
  }
}
