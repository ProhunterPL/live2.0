"""
LIVE 2.0 Cluster Watcher - Monitor clustershots folder for new images.

This script watches the clustershots/ directory and automatically processes
new cluster images by calling matcher.py.

Usage:
    python matcher/watcher.py

The watcher will:
1. Monitor clustershots/ for new PNG/JPG files
2. Wait for corresponding JSON file to appear
3. Run matcher.py to generate comparison panel
4. Continue monitoring for new files
"""
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time
from pathlib import Path
import subprocess
import sys


WATCH_DIR = "clustershots"
EXTENSIONS = {".png", ".jpg", ".jpeg"}


class ClusterHandler(FileSystemEventHandler):
    """File system event handler for cluster images."""
    
    def __init__(self):
        super().__init__()
        self.processing = set()
    
    def on_created(self, event):
        """Handle file creation events."""
        if event.is_directory:
            return
        
        p = Path(event.src_path)
        
        # Only process image files
        if p.suffix.lower() not in EXTENSIONS:
            return
        
        # Avoid processing the same file multiple times
        if str(p) in self.processing:
            return
        
        self.processing.add(str(p))
        
        print(f"\n🔔 New cluster image detected: {p.name}")
        
        # Wait for corresponding JSON file
        json_path = p.with_suffix(".json")
        print(f"⏳ Waiting for metadata: {json_path.name}")
        
        max_wait = 40  # 4 seconds max wait
        for i in range(max_wait):
            if json_path.exists():
                print(f"✓ Found metadata after {i*0.1:.1f}s")
                break
            time.sleep(0.1)
        else:
            print(f"⚠️  Warning: JSON metadata not found after {max_wait*0.1}s")
            print(f"   Skipping {p.name}")
            self.processing.discard(str(p))
            return
        
        # Process the pair
        print(f"🚀 Starting matcher for {p.name}...")
        try:
            result = subprocess.run(
                [sys.executable, "matcher/matcher.py", str(p)],
                check=True,
                capture_output=True,
                text=True
            )
            print(result.stdout)
            print("✅ Match completed successfully!")
        except subprocess.CalledProcessError as e:
            print(f"❌ Matcher failed with exit code {e.returncode}")
            print("STDOUT:", e.stdout)
            print("STDERR:", e.stderr)
        except Exception as e:
            print(f"❌ Unexpected error running matcher: {e}")
        finally:
            self.processing.discard(str(p))


def main():
    """Start the watcher."""
    # Ensure watch directory exists
    watch_path = Path(WATCH_DIR)
    watch_path.mkdir(exist_ok=True)
    
    print("="*60)
    print("LIVE 2.0 Cluster Watcher")
    print("="*60)
    print(f"📁 Watching: {watch_path.absolute()}")
    print(f"🔍 Looking for: {', '.join(EXTENSIONS)}")
    print(f"💡 Save cluster images + JSON to this folder")
    print(f"⚙️  Press Ctrl+C to stop")
    print("="*60)
    
    # Create observer
    observer = Observer()
    handler = ClusterHandler()
    observer.schedule(handler, str(watch_path), recursive=False)
    observer.start()
    
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("\n\n🛑 Stopping watcher...")
        observer.stop()
        observer.join()
        print("✓ Watcher stopped cleanly")


if __name__ == "__main__":
    main()

