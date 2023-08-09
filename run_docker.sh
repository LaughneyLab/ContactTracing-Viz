#!/usr/bin/env sh

# Add -d to daemonize
( docker stop ct_viz || true ) && ( docker rm ct_viz || true )
docker run --rm --name ct_viz -p 8000:8000 -v data:/app/data ct-site