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
{% end if %}
{% end for %}

---

## Geometry of Polynomials
{% for pub in site.publications reversed %}
{% if pub.tags contains "Geometry of Polynomials" %}
  {% include archive-single.html type="list" %}
{% end if %}
{% end for %}

---

## Matrix Inequalities
{% for pub in site.publications reversed %}
{% if pub.tags contains "Matrix Inequalities" %}
  {% include archive-single.html type="list" %}
{% end if %}
{% end for %}
