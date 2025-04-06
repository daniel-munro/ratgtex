set -e

bundle exec jekyll build
cd _site
ln -s ../data data
echo "Run 'cd _site/api && python3 app.py' in another tab to start the flask API server (for testing only)."
# Redirect /api/ calls to flask test server:
http-server -p 4000 -P http://127.0.0.1:5000
