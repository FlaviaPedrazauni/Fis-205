# 🍞 Simulador Termodinámico de Horneado de Pan

**Efecto de las condiciones de horneado en la transferencia de calor**

![Python](https://img.shields.io/badge/Python-3.13-blue?logo=python&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-013243?logo=numpy&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-8CAAE6?logo=scipy&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-11557c)
![ipywidgets](https://img.shields.io/badge/ipywidgets-interactive-orange)
![Estado](https://img.shields.io/badge/estado-funcional-brightgreen)

Proyecto del curso **Fis-205** (Universidad Técnica Federico Santa María). Se modela la
transferencia de calor en una masa de pan durante el horneado mediante un **solver 1D de la
ecuación del calor** resuelto con el método de **Crank–Nicolson**. El objetivo es estudiar cómo
las condiciones del horno afectan el perfil térmico interno, la formación de corteza y la
pérdida de humedad, comparando los resultados cualitativamente con la literatura (Zürcher y
Purlis).

El proyecto incluye, además del solver, una **animación (GIF)** del perfil térmico evolucionando
en el tiempo y un **dashboard interactivo** (sliders con `ipywidgets`) para explorar el efecto de
los parámetros de horneado en tiempo real.

---

## 📑 Tabla de contenidos

- [Descripción del modelo](#-descripción-del-modelo)
  - [Fase 1 — Propiedades constantes](#fase-1--propiedades-constantes-geometría-plana)
  - [Fase 2 — Propiedades variables + evaporación](#fase-2--propiedades-variables--evaporación-geometría-cilíndrica)
- [Visualizaciones](#-visualizaciones-de-propiedades-térmicas)
  - [Animación del perfil térmico (GIF)](#-animación-del-perfil-térmico-gif)
  - [Dashboard interactivo](#-dashboard-interactivo)
- [Estructura del repositorio](#-estructura-del-repositorio)
- [Requisitos e instalación](#-requisitos-e-instalación)
- [Uso](#-uso)
- [Resultados](#-resultados)
- [Referencias](#-referencias)
- [Autora](#-autora)

---

## 🔬 Descripción del modelo

El pan se trata como un **medio sólido húmedo** en el que la difusión de calor domina el proceso.
La simulación está dividida en dos fases de complejidad creciente.

### Fase 1 — Propiedades constantes (geometría plana)

Primer modelo de validación: ecuación del calor 1D en coordenadas cartesianas con propiedades
termofísicas **constantes**.

$$\rho\, C_p\, \frac{\partial T}{\partial t} = k\, \frac{\partial^2 T}{\partial x^2}, \qquad \alpha = \frac{k}{\rho\, C_p}$$

**Condiciones de frontera:**

- **Centro ($x = 0$):** simetría / adiabática, $\left.\dfrac{\partial T}{\partial x}\right|_{x=0} = 0$
- **Superficie ($x = L$):** condición de **Robin** (convección + radiación)

$$-k\,\left.\frac{\partial T}{\partial x}\right|_{x=L} = H\,\bigl(T_{\text{sup}} - T_{\text{horno}}\bigr)$$

La discretización temporal usa Crank–Nicolson con el coeficiente $r = \dfrac{\alpha\,\Delta t}{2\,\Delta x^2}$,
generando un **sistema tridiagonal** que se resuelve eficientemente con `scipy.linalg.solve_banded`.

| Parámetro | Símbolo | Valor |
|---|---|---|
| Espesor (centro → superficie) | $L$ | 0.04 m |
| Nodos espaciales | $N_x$ | 140 |
| Paso temporal | $\Delta t$ | 0.5 s |
| Tiempo de horneado | $t_{\max}$ | 1800 s (30 min) |
| Densidad | $\rho$ | 800 kg/m³ |
| Calor específico | $C_p$ | 2000 J/(kg·K) |
| Conductividad térmica | $k$ | 0.5 W/(m·K) |
| Temperatura del horno | $T_{\text{horno}}$ | 180 °C |
| Coef. convección + radiación | $H$ | 25 W/(m²·K) |
| Temperatura inicial | $T_{\text{ini}}$ | 25 °C |

### Fase 2 — Propiedades variables + evaporación (geometría cilíndrica)

El pan se modela como un **cilindro infinito** con propiedades que **dependen de la temperatura y
la humedad**, e incorpora un **término sumidero $S$** que representa el calor latente consumido por
la evaporación del agua interna.

$$\rho(T)\, C_p(T, W)\, \frac{\partial T}{\partial t} = \frac{1}{r}\,\frac{\partial}{\partial r}\!\left(r\, k(T)\, \frac{\partial T}{\partial r}\right) - S$$

Aspectos numéricos y físicos implementados:

- **Propiedades efectivas** `k_eff(T)`, `Cp_sens(T, W)` , `H_eff(Ts)`y `rho_eff(T)`, tomadas de Purlis, que
  cambian entre las ramas de **miga** y **corteza**.
- **Heaviside suavizada** `smoothH(T)`: transición continua (≈0 bajo 100 °C, ≈1 sobre 100 °C)
  centrada en $T_f$, que evita saltos numéricos al mezclar ambas ramas.
- **Cambio de fase por *operator-splitting*** (`factor_evaporacion`): la energía por encima de
  100 °C se invierte en evaporar agua (calor latente $L_v = 2.3\times10^6$ J/kg), manteniendo el
  nodo en la **meseta de 100 °C** hasta que se seca y vuelve a recalentarse.
- **Tratamiento del centro** mediante L'Hôpital (factor 2) para evitar la singularidad $1/r$ en $r = 0$.
- **Sub-pasos temporales** (`n_sub = 3`) para estabilizar las oscilaciones del *operator-splitting*
  y suavizar las curvas de temperatura.

Constantes del cambio de fase: $T_f = 373.16$ K (100 °C), semiancho de banda $\Delta T_b = 0.5$ K,
densidad del sólido seco $\rho_s = 241.76$ kg/m³, humedad inicial $W_{\text{ini}} = 0.60$ kg agua/kg sólido.

---

## 🎬 Visualizaciones interactivas

Además de los gráficos estáticos, el notebook genera dos productos de visualización para explorar
la dinámica del horneado.

### 🎞️ Animación del perfil térmico (GIF)

Se genera una **animación** que muestra cómo evoluciona el perfil de temperatura $T(x)$ a lo largo
de todo el horneado: la superficie supera rápidamente los 100 °C (frente de evaporación y formación
de corteza) mientras el centro se calienta de forma progresiva y desarrolla la meseta característica.

La animación se construye con `matplotlib.animation.FuncAnimation` y se exporta a GIF con el
escritor `pillow`, recorriendo los pasos guardados del bucle temporal. El archivo resultante se
guarda en `figuras/animacion_horneado.gif`.

```python
from matplotlib import animation
# ... tras correr la simulación y guardar los perfiles T(x) por paso ...
ani = animation.FuncAnimation(fig, update, frames=n_frames, interval=40, blit=True)
ani.save("figuras/animacion_horneado.gif", writer="pillow", fps=25)
```

### 🎛️ Dashboard interactivo

Un **panel interactivo con `ipywidgets`** permite variar los parámetros de horneado mediante
*sliders* y volver a graficar las curvas $T(t)$ del centro y la superficie sin editar el código.
Parámetros expuestos:

- 🌡️ Temperatura del horno $T_{\text{horno}}$
- 💨 Coeficiente convección 
- 💧 Humedad inicial $W_{\text{ini}}$
- 📏 Espesor de la pieza $L$

El dashboard usa `interact` / `interactive` de `ipywidgets`, de modo que cada cambio en un slider
re-ejecuta el solver y actualiza los gráficos en vivo dentro del notebook.

> ℹ️ El dashboard es interactivo solo al **ejecutar el notebook** (Jupyter / JupyterLab). En la
> vista previa de GitHub los widgets aparecen estáticos; para usarlos hay que correr las celdas
> localmente o en un entorno como Binder/Colab con `ipywidgets` habilitado.

---

## 📂 Estructura del repositorio

```text
Fis-205/
├── README.md
├── requirements.txt
├── notebooks/
│   └── simulador_horneado.ipynb      # código completo (Fase 1 + Fase 2 + GIF)
|   └── Dashboard_horneado(1d).ipynb  # código completo (Panel interactivo)
├── informe/
│   └── informe_avance.pdf            # informe del proyecto
└── figuras/
    ├── animacion_horneado.gif        # animación del perfil térmico
    ├── perfil_temperatura.png        # gráficas estáticas exportadas
    └── ...
```

---

## ⚙️ Requisitos e instalación

El proyecto usa Python 3 y librerías científicas estándar. El GIF requiere `pillow` y el panel
interactivo requiere `ipywidgets`.

```bash
# Clonar el repositorio
git clone https://github.com/FlaviaPedrazauni/Fis-205.git
cd Fis-205

# (Opcional) crear un entorno virtual
python -m venv .venv
source .venv/bin/activate      # En Windows: .venv\Scripts\activate

# Instalar dependencias
pip install -r requirements.txt
```

**`requirements.txt`:**

```text
numpy
scipy
matplotlib
ipywidgets
pillow
jupyter
```

---

## ▶️ Uso

Abre el notebook con Jupyter y ejecuta las celdas en orden:

```bash
jupyter notebook notebooks/simulador_horneado.ipynb
```

El notebook está organizado de forma secuencial:
-**Código principal:**
1. **Parámetros y malla** — geometría, tiempo y propiedades termofísicas.
2. **Fase 1** — construcción de la matriz tridiagonal, bucle temporal y visualización del perfil
   espacial y la evolución temporal (centro vs. superficie).
3. **Fase 2** — definición de propiedades variables, ensamblaje del sistema cilíndrico, término de
   evaporación, bucle principal y los tres gráficos de resultados.
4. **Animación (GIF)** — genera y exporta `figuras/animacion_horneado.gif`.
-**Dashboard interactivo:** Sliders de `ipywidgets` para variar $T_{\text{horno}}$, $H$,
   $W_{\text{ini}}$ y $L$, y observar las curvas $T(t)$ en tiempo real.

> Para que los sliders funcionen, ejecuta el notebook dashboard_horneado(1d) (no basta con la vista previa de
> GitHub) y asegúrate de tener `ipywidgets` instalado y habilitado.

---

## 📊 Resultados

**Fase 1.** La superficie se calienta rápidamente y supera los 100 °C (inicio de la formación de
corteza), mientras que el centro alcanza solo ~59 °C tras 30 minutos. Esto coincide
cualitativamente con el modelo de propiedades constantes de Zürcher y motiva el paso a la Fase 2:
para que el centro llegue a ~99 °C se requiere extender el horneado a ~60 min.

**Fase 2.** Al incorporar propiedades variables y el calor latente de evaporación, el centro
desarrolla la **meseta característica a 100 °C** (cocción y evaporación), la **corteza se forma en
la superficie** y el modelo entrega tres salidas:

- 🌡️ **Temperatura** del centro y la superficie en función del tiempo.
- ⚖️ **Pérdida de peso** (%) en función del tiempo de horneado.
- 💧 **Humedad final** en función del radio (frente de secado hacia la superficie).

**Visualizaciones.** La **animación** permite ver el avance del frente térmico y la meseta de
evaporación cuadro a cuadro, mientras que el **dashboard** muestra de forma inmediata cómo un horno
más caliente, un mayor $H$ o un menor espesor aceleran el calentamiento del centro y la formación
de corteza.

---

## 📚 Referencias

1. Purlis, E., & Salvadori, V. O. (2009a). *Bread baking as a moving boundary problem. Part 1: Model development.* Journal of Food Engineering, 91(3), 428–433.
2. Purlis, E., & Salvadori, V. O. (2009b). *Bread baking as a moving boundary problem. Part 2: Model validation and numerical simulation.* Journal of Food Engineering, 91(3), 434–441.
3. Purlis, E. (2011). *Bread baking: Technological considerations based on process modelling and simulation.* Journal of Food Engineering, 103(1), 92-102.
4. Sluimer, P., & Krist-Spit, C. E. (1987). *Some aspects of heat and mass transfer in baking.* Proceedings of the 11th International Congress of Mechanical Engineering, 1, 55–64.
5. Nicolas, V., Glouannec, P., Ploteau, J. P., Jury, V., & Le-Bail, A. (2016). *A multi-physics model for bread baking including heat and mass transfer, CO₂ production, and volume expansion.* International Journal of Heat and Mass Transfer, 98, 223–235.
6. Zheleva, I., & Kambourova, V. (2005). *Identification of heat and mass transfer processes in bread during baking.* Thermal Science, 9(2), 73–86.
7. Zürcher, S. (2014). *A two-state model for bread baking.* Food and Bioproducts Processing, 92(1), 16–23.
8. Hamdami, N., Monteau, J.-Y., and Le Bail, A. (2004).*Thermophysical properties evolution of French partly baked bread during freezing.* Food Research International, 37(7):703--713.

---

## ✍️ Autora

**Flavia Pedraza Alarcón** — Universidad Técnica Federico Santa María
📧 fpedraza@usm.cl

*Proyecto académico del curso Fis-205.*
