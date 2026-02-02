/**
 * Theme management for EVN Observation Planner
 * Handles dark/light mode switching with system preference detection
 */

(function() {
    'use strict';

    const THEME_KEY = 'planobs-theme';

    /**
     * Get the current effective theme
     * @returns {'light' | 'dark'}
     */
    function getEffectiveTheme() {
        const stored = localStorage.getItem(THEME_KEY);
        if (stored === 'light' || stored === 'dark') {
            return stored;
        }
        // Check system preference
        if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
            return 'dark';
        }
        return 'light';
    }

    /**
     * Apply theme to document
     * @param {'light' | 'dark'} theme
     */
    function applyTheme(theme) {
        document.documentElement.setAttribute('data-theme', theme);
        
        // Update Plotly charts if they exist
        updatePlotlyTheme(theme);
        
        // Update Mantine provider if available
        updateMantineTheme(theme);
    }

    /**
     * Update Plotly chart colors for theme
     * @param {'light' | 'dark'} theme
     */
    function updatePlotlyTheme(theme) {
        const isDark = theme === 'dark';
        const textColor = isDark ? '#e9ecef' : '#344767';
        const gridColor = isDark ? '#495057' : '#e0e0e0';
        const bgColor = 'rgba(0,0,0,0)';

        // Find all Plotly charts and update their layout
        const plots = document.querySelectorAll('.js-plotly-plot');
        plots.forEach(function(plot) {
            if (plot._fullLayout) {
                try {
                    Plotly.relayout(plot, {
                        'paper_bgcolor': bgColor,
                        'plot_bgcolor': bgColor,
                        'font.color': textColor,
                        'title.font.color': textColor,
                        'xaxis.gridcolor': gridColor,
                        'xaxis.linecolor': textColor,
                        'xaxis.tickcolor': textColor,
                        'xaxis.tickfont.color': textColor,
                        'xaxis.title.font.color': textColor,
                        'xaxis2.linecolor': textColor,
                        'xaxis2.tickcolor': textColor,
                        'xaxis2.tickfont.color': textColor,
                        'xaxis2.title.font.color': textColor,
                        'yaxis.gridcolor': gridColor,
                        'yaxis.linecolor': textColor,
                        'yaxis.tickcolor': textColor,
                        'yaxis.tickfont.color': textColor,
                        'yaxis.title.font.color': textColor,
                        'legend.font.color': textColor,
                        'coloraxis.colorbar.tickfont.color': textColor,
                        'coloraxis.colorbar.title.font.color': textColor,
                        'geo.bgcolor': bgColor
                    });
                } catch (e) {
                    // Plotly may not be fully initialized yet
                    console.debug('Could not update plot theme:', e);
                }
            }
        });
    }

    /**
     * Update Mantine theme if MantineProvider is available
     * @param {'light' | 'dark'} theme
     */
    function updateMantineTheme(theme) {
        // Mantine components will pick up CSS variables automatically
        // but we can add a class for additional styling hooks
        document.body.classList.remove('theme-light', 'theme-dark');
        document.body.classList.add('theme-' + theme);
    }

    /**
     * Toggle between light and dark themes
     * @returns {'light' | 'dark'} The new theme
     */
    function toggleTheme() {
        const current = getEffectiveTheme();
        const newTheme = current === 'dark' ? 'light' : 'dark';
        localStorage.setItem(THEME_KEY, newTheme);
        applyTheme(newTheme);
        return newTheme;
    }

    /**
     * Initialize theme on page load
     */
    function initTheme() {
        const theme = getEffectiveTheme();
        applyTheme(theme);

        // Listen for system theme changes
        if (window.matchMedia) {
            window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', function(e) {
                // Only apply if user hasn't manually set a preference
                if (!localStorage.getItem(THEME_KEY)) {
                    applyTheme(e.matches ? 'dark' : 'light');
                }
            });
        }
    }

    // Expose functions globally for Dash callbacks
    window.planobsTheme = {
        toggle: toggleTheme,
        get: getEffectiveTheme,
        set: function(theme) {
            localStorage.setItem(THEME_KEY, theme);
            applyTheme(theme);
        },
        init: initTheme,
        updatePlots: function() {
            updatePlotlyTheme(getEffectiveTheme());
        }
    };

    // Initialize on DOM ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', initTheme);
    } else {
        initTheme();
    }

    // Also update plots when they're created/updated
    // MutationObserver to watch for new plots and plot updates
    let plotUpdateTimeout = null;
    const observer = new MutationObserver(function(mutations) {
        let shouldUpdate = false;
        mutations.forEach(function(mutation) {
            // Check for new plot nodes
            mutation.addedNodes.forEach(function(node) {
                if (node.nodeType === 1) {  // Element node
                    if (node.classList && node.classList.contains('js-plotly-plot')) {
                        shouldUpdate = true;
                    }
                    // Also check children for plots
                    if (node.querySelector && node.querySelector('.js-plotly-plot')) {
                        shouldUpdate = true;
                    }
                }
            });
            // Check for modifications to existing plots
            if (mutation.target && mutation.target.closest && mutation.target.closest('.js-plotly-plot')) {
                shouldUpdate = true;
            }
        });
        
        if (shouldUpdate) {
            // Debounce plot updates
            if (plotUpdateTimeout) clearTimeout(plotUpdateTimeout);
            plotUpdateTimeout = setTimeout(function() {
                updatePlotlyTheme(getEffectiveTheme());
            }, 200);
        }
    });

    observer.observe(document.body, { childList: true, subtree: true, attributes: true });

    // Also hook into Dash's callback completion if available
    if (window.dash_clientside) {
        const originalNoUpdate = window.dash_clientside.no_update;
        // Periodically check for new plots after any Dash activity
        setInterval(function() {
            const plots = document.querySelectorAll('.js-plotly-plot');
            if (plots.length > 0) {
                updatePlotlyTheme(getEffectiveTheme());
            }
        }, 2000);
    }

})();
