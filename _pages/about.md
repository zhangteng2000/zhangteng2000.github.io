---
title: "About me"
permalink: /
excerpt: "About me"
author_profile: true
redirect_from:
  - /about/
  - /about.html
---

I am a Ph.D. student in Mathematics at the  [School of Mathematics and Statistics](https://math.xjtu.edu.cn/),  [Xi’an Jiaotong University, China](https://www.xjtu.edu.cn/),  under the supervision of  [Prof. Minghua Lin](http://gr.xjtu.edu.cn/en/web/mh.lin) and a joint Ph.D. student at the [Department of Analysis, Bolyai Institute](https://www.math.u-szeged.hu/analysis/),  [University of Szeged, Hungary](https://u-szeged.hu/),  under the supervision of  [Prof. Lajos Molnár](https://www.math.u-szeged.hu/~molnarl/).


My research interests lie in functional analysis, operator theory,  and geometry of polynomials. More specifically, I study  *transformations on operator algebras*,  *function algebras*, and *quantum structures*,  with a focus on *preserver problems*,  *local and nonlinear transformations*,  *matrix inequalities* and *value distribution problems of polynomials*.


📧 Email: [teng.zhang@stu.xjtu.edu.cn](mailto:teng.zhang@stu.xjtu.edu.cn)


## Education

- **Ph.D. in Mathematics**,  
  Xi'an Jiaotong University, Xi'an, China, *09/2021–Present*  
  School of Mathematics and Statistics  
  Advisor: Minghua Lin

- **Joint Ph.D. Student in Mathematics**,  
  University of Szeged, Szeged, Hungary, *12/2025–12/2026*  
  Department of Analysis, Bolyai Institute  
  Advisor: Lajos Molnár


- **B.S. in Mathematics**,  
  Taiyuan University of Technology, Taiyuan, China, *09/2017–07/2021*

---

## Publications

{% raw %}{% include base_path %}

### Preserver Problems
{% for pub in site.publications reversed %}
{% if pub.tags contains "Preserver" %}
  {% include archive-single.html type="list" %}
{% endif %}
{% endfor %}

---

### Geometry of Polynomial
{% for pub in site.publications reversed %}
{% if pub.tags contains "Geometry of Polynomials" %}
  {% include archive-single.html type="list" %}
{% endif %}
{% endfor %}

---

### Matrix Inequalities
{% for pub in site.publications reversed %}
{% if pub.tags contains "Matrix Inequalities" %}
  {% include archive-single.html type="list" %}
{% endif %}
{% endfor %}
{% endraw %}


