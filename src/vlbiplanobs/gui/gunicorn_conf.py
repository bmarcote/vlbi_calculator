"""Gunicorn settings to serve the EVN Observation Planner web app in production.

Usage (works from any directory, the settings are loaded from the installed package):

    gunicorn -c python:vlbiplanobs.gui.gunicorn_conf vlbiplanobs.gui.main:server

Every value can be changed through the environment variables named below, without editing this file.
"""
import os

bind = os.environ.get('PLANOBS_BIND', '127.0.0.1:8050')

# The app is imported (and warmed up, see vlbiplanobs.gui.main.warm_up) once in the master process and the
# workers are forked from it. Hence all workers answer their very first request at full speed, a worker
# that is replaced (max_requests, timeout) is ready immediately, and the memory of the loaded libraries and
# catalogs is shared between the workers instead of being duplicated in each of them.
preload_app = True

# Threaded workers: a worker computing an observation or waiting for a PDF to be rendered can still serve
# the (many, tiny) requests of static files of other users. numpy/erfa release the GIL in the heavy parts.
worker_class = 'gthread'
workers = int(os.environ.get('PLANOBS_WORKERS', min(4, os.cpu_count() or 1)))
threads = int(os.environ.get('PLANOBS_THREADS', 4))

# PDF reports with figures can take several seconds.
timeout = int(os.environ.get('PLANOBS_TIMEOUT', 120))
graceful_timeout = 30
keepalive = 5

# Replace the workers every now and then to keep the memory bounded (cheap, thanks to preload_app).
max_requests = int(os.environ.get('PLANOBS_MAX_REQUESTS', 5000))
max_requests_jitter = max_requests // 10

# Each worker keeps its heartbeat file in memory instead of on disk.
worker_tmp_dir = '/dev/shm' if os.path.isdir('/dev/shm') else None
