import '@testing-library/jest-dom'

// Mock WebSocket
global.WebSocket = class WebSocket {
  constructor(public url: string) {}
  
  close() {}
  send() {}
  
  onopen = null
  onclose = null
  onmessage = null
  onerror = null
  
  readyState = WebSocket.OPEN
} as any

// Mock msgpack
jest.mock('msgpack-lite', () => ({
  encode: jest.fn((data) => JSON.stringify(data)),
  decode: jest.fn((data) => JSON.parse(data))
}))

// Mock canvas
HTMLCanvasElement.prototype.getContext = jest.fn(() => ({
  fillRect: jest.fn(),
  fillStyle: '',
  strokeStyle: '',
  lineWidth: 0,
  beginPath: jest.fn(),
  moveTo: jest.fn(),
  lineTo: jest.fn(),
  stroke: jest.fn(),
  arc: jest.fn(),
  fill: jest.fn(),
  scale: jest.fn(),
  getBoundingClientRect: jest.fn(() => ({
    width: 800,
    height: 600,
    left: 0,
    top: 0
  }))
}))

// Mock ResizeObserver
global.ResizeObserver = class ResizeObserver {
  constructor(callback: ResizeObserverCallback) {}
  observe() {}
  unobserve() {}
  disconnect() {}
}
