---
title: "Publications"
permalink: /publications/
layout: single
author_profile: true
---

{% include base_path %}

## Preserver Problems
{% for pub in site.publications reversed %}
{% if pub.tags contains "Preserver" %}
  {% include archive-single.html type="list" %}
{% endif %}
{% endfor %}

---

## Polynomial Geometry
{% for pub in site.publications reversed %}
{% if pub.tags contains "Polynomial Geometry" %}
  {% include archive-single.html type="list" %}
{% endif %}
{% endfor %}

---

## Norm Inequalities
{% for pub in site.publications reversed %}
{% if pub.tags contains "Norm Inequalities" %}
  {% include archive-single.html type="list" %}
{% endif %}
{% endfor %}
