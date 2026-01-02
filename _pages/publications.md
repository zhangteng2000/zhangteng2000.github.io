---
title: "Publications"
permalink: /publications/
layout: single
author_profile: true
---

### DEBUG: All publications

{% for pub in site.publications %}
- {{ pub.title }} | tags={{ pub.tags }} | year={{ pub.year }}
{% endfor %}




