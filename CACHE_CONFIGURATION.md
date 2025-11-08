# Cache Configuration for Optimal Performance

## Overview
This document provides cache configuration recommendations for hosting the website on various platforms.

## GitHub Pages (Current Platform)

GitHub Pages automatically sets cache headers, but you cannot customize them. The default headers are:
- **HTML files**: `Cache-Control: max-age=600` (10 minutes)
- **Assets (CSS, JS, images)**: `Cache-Control: max-age=3600` (1 hour)

### Workaround for Better Caching
To achieve better caching on GitHub Pages, use versioned filenames:
- `enhanced-styles.min.css?v=1.0.0`
- `enhanced-scripts.min.js?v=1.0.0`

When you update files, increment the version number to bust the cache.

## Netlify Configuration

If migrating to Netlify, create `netlify.toml` in the root:

```toml
[[headers]]
  for = "/*.html"
  [headers.values]
    Cache-Control = "public, max-age=0, must-revalidate"
    X-Content-Type-Options = "nosniff"
    X-Frame-Options = "DENY"
    X-XSS-Protection = "1; mode=block"
    Referrer-Policy = "strict-origin-when-cross-origin"

[[headers]]
  for = "/*.css"
  [headers.values]
    Cache-Control = "public, max-age=31536000, immutable"

[[headers]]
  for = "/*.js"
  [headers.values]
    Cache-Control = "public, max-age=31536000, immutable"

[[headers]]
  for = "/images/*"
  [headers.values]
    Cache-Control = "public, max-age=31536000, immutable"

[[headers]]
  for = "/*.woff2"
  [headers.values]
    Cache-Control = "public, max-age=31536000, immutable"
```

## Vercel Configuration

If migrating to Vercel, create `vercel.json` in the root:

```json
{
  "headers": [
    {
      "source": "/(.*).html",
      "headers": [
        {
          "key": "Cache-Control",
          "value": "public, max-age=0, must-revalidate"
        }
      ]
    },
    {
      "source": "/(.*).css",
      "headers": [
        {
          "key": "Cache-Control",
          "value": "public, max-age=31536000, immutable"
        }
      ]
    },
    {
      "source": "/(.*).js",
      "headers": [
        {
          "key": "Cache-Control",
          "value": "public, max-age=31536000, immutable"
        }
      ]
    },
    {
      "source": "/images/(.*)",
      "headers": [
        {
          "key": "Cache-Control",
          "value": "public, max-age=31536000, immutable"
        }
      ]
    }
  ]
}
```

## Apache (.htaccess)

If using Apache server, create/update `.htaccess`:

```apache
<IfModule mod_expires.c>
    ExpiresActive On

    # HTML files
    ExpiresByType text/html "access plus 0 seconds"

    # CSS and JavaScript
    ExpiresByType text/css "access plus 1 year"
    ExpiresByType application/javascript "access plus 1 year"
    ExpiresByType application/x-javascript "access plus 1 year"

    # Images
    ExpiresByType image/jpeg "access plus 1 year"
    ExpiresByType image/png "access plus 1 year"
    ExpiresByType image/gif "access plus 1 year"
    ExpiresByType image/svg+xml "access plus 1 year"
    ExpiresByType image/webp "access plus 1 year"

    # Fonts
    ExpiresByType font/woff2 "access plus 1 year"
    ExpiresByType application/font-woff2 "access plus 1 year"
</IfModule>

<IfModule mod_headers.c>
    # Security headers
    Header set X-Content-Type-Options "nosniff"
    Header set X-Frame-Options "DENY"
    Header set X-XSS-Protection "1; mode=block"
    Header set Referrer-Policy "strict-origin-when-cross-origin"

    # Cache control for static assets with versioning
    <FilesMatch "\.(css|js|jpg|jpeg|png|gif|svg|webp|woff2)$">
        Header set Cache-Control "public, max-age=31536000, immutable"
    </FilesMatch>

    # No cache for HTML
    <FilesMatch "\.(html|htm)$">
        Header set Cache-Control "public, max-age=0, must-revalidate"
    </FilesMatch>
</IfModule>
```

## Nginx Configuration

If using Nginx, add to your server block:

```nginx
# Cache static assets for 1 year
location ~* \.(css|js|jpg|jpeg|png|gif|svg|webp|woff2)$ {
    expires 1y;
    add_header Cache-Control "public, immutable";
}

# No cache for HTML
location ~* \.(html|htm)$ {
    expires off;
    add_header Cache-Control "public, max-age=0, must-revalidate";
}

# Security headers
add_header X-Content-Type-Options "nosniff" always;
add_header X-Frame-Options "DENY" always;
add_header X-XSS-Protection "1; mode=block" always;
add_header Referrer-Policy "strict-origin-when-cross-origin" always;
```

## Cloudflare CDN Setup

If using Cloudflare as a CDN:

1. **Enable Cloudflare Caching**:
   - Go to Cloudflare Dashboard → Caching → Configuration
   - Set Browser Cache TTL: 1 year
   - Enable "Respect Existing Headers"

2. **Page Rules** (create these):
   - Pattern: `*ajaykhanna.github.io/*.css*`
     - Cache Level: Cache Everything
     - Edge Cache TTL: 1 month
     - Browser Cache TTL: 1 year

   - Pattern: `*ajaykhanna.github.io/*.js*`
     - Cache Level: Cache Everything
     - Edge Cache TTL: 1 month
     - Browser Cache TTL: 1 year

   - Pattern: `*ajaykhanna.github.io/images/*`
     - Cache Level: Cache Everything
     - Edge Cache TTL: 1 year
     - Browser Cache TTL: 1 year

3. **Enable Auto Minify**:
   - Go to Speed → Optimization
   - Enable: JavaScript, CSS, HTML

4. **Enable Brotli Compression**:
   - Go to Speed → Optimization
   - Enable Brotli compression

## Recommended Cache Strategy

### Cache Duration by File Type:

| File Type | Cache Duration | Reasoning |
|-----------|----------------|-----------|
| HTML | No cache / 0s | Content changes frequently |
| CSS/JS (minified) | 1 year | Immutable, versioned filenames |
| Images | 1 year | Rarely change |
| Fonts | 1 year | Never change |
| JSON data | 1 hour | May update periodically |
| API responses | Variable | Depends on data volatility |

### Cache-Control Header Values:

- **`public`**: Can be cached by browsers and CDNs
- **`private`**: Only browser cache, not CDNs
- **`max-age=31536000`**: Cache for 1 year (in seconds)
- **`immutable`**: File will never change (perfect for versioned assets)
- **`must-revalidate`**: Check with server before using stale cache
- **`no-cache`**: Always revalidate with server
- **`no-store`**: Never cache

## Version-Based Cache Busting

To ensure users get the latest files when you update, use version parameters:

### In HTML:
```html
<link rel="stylesheet" href="enhanced-styles.min.css?v=2.0.0">
<script src="enhanced-scripts.min.js?v=2.0.0"></script>
```

### Update Process:
1. Make changes to CSS/JS
2. Increment version number in HTML
3. Commit and push
4. Users automatically get new version

### Automated Versioning (Optional):

Use git commit hash as version:
```bash
# Generate version from git
VERSION=$(git rev-parse --short HEAD)
sed -i "s/\?v=[^\"]*/?v=$VERSION/g" index.html
```

## Performance Testing

After implementing caching, test with:

1. **Google PageSpeed Insights**: https://pagespeed.web.dev/
2. **GTmetrix**: https://gtmetrix.com/
3. **WebPageTest**: https://www.webpagetest.org/
4. **Chrome DevTools**:
   - Network tab → Check cache status
   - Lighthouse → Performance audit

## Current Setup (GitHub Pages)

The website is currently hosted on GitHub Pages which:
- ✅ Automatically serves content via CDN
- ✅ Has basic caching (1 hour for assets)
- ✅ Supports HTTPS
- ❌ Doesn't allow custom cache headers
- ❌ Doesn't support edge functions

## Future Optimization Opportunities

If you migrate from GitHub Pages, consider:

1. **Cloudflare** (Free tier):
   - Global CDN
   - Custom cache rules
   - Auto minification
   - Brotli compression

2. **Netlify** (Free tier):
   - Custom headers
   - Asset optimization
   - Instant cache invalidation
   - Branch previews

3. **Vercel** (Free tier):
   - Edge caching
   - Auto-optimization
   - Analytics
   - Fast builds

## Monitoring Cache Performance

Track cache effectiveness with:

1. **Cache Hit Ratio**: Aim for >90%
2. **Bandwidth Savings**: Monitor reduced origin requests
3. **Load Times**: Compare cached vs uncached loads
4. **User Experience**: Core Web Vitals metrics

---

**Last Updated**: November 8, 2025
**Status**: Optimized for GitHub Pages with versioning strategy
