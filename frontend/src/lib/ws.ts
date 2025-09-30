// Live 2.0 WebSocket Client

import type { SimulationData } from './types'
// eslint-disable-next-line @typescript-eslint/ban-ts-comment
// @ts-ignore - no types for msgpack-lite
import msgpack from 'msgpack-lite'

export class WebSocketClient {
  private ws: WebSocket | null = null
  private url: string
  private reconnectAttempts = 0
  private maxReconnectAttempts = 5
  private reconnectDelay = 1000
  private listeners: Map<string, Function[]> = new Map()

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
          
          // Only reconnect if it's not a normal closure and we have a valid simulation ID
          if (event.code !== 1000 && this.reconnectAttempts < this.maxReconnectAttempts && urlOrId) {
            this.reconnect(urlOrId)
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
    if (this.ws) {
      this.ws.close(1000, 'Client disconnect')
      this.ws = null
    }
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
