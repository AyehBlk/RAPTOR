"""
RAPTOR Diagnostic Suite — Comprehensive health check and debugging tool.

Run routinely before releases, after adding new features/modules, or when
something breaks. Covers structure, imports, dependencies, tests, CLI,
dashboard, code quality, git status, and more.

Usage:
    python check_raptor.py                  # Full diagnostic (all 16 checks)
    python check_raptor.py --quick          # Skip slow checks (tests, CLI)
    python check_raptor.py --fix            # Show fix suggestions for every failure
    python check_raptor.py --section 3      # Run only section 3 (imports)
    python check_raptor.py --section 3,5,7  # Run sections 3, 5, and 7
    python check_raptor.py --verbose        # Show extra debug info
    python check_raptor.py --save           # Save results to check_results.txt

How to add a new module:
    1. Add the package path to EXPECTED_PACKAGES
    2. Add module paths to EXPECTED_CORE_MODULES or EXPECTED_ACQUISITION_MODULES
    3. Add expected exports to EXPECTED_RAPTOR_EXPORTS or EXPECTED_ACQUISITION_EXPORTS
    4. Add any new required/optional deps to REQUIRED_DEPS or OPTIONAL_DEPS
    5. Update EXPECTED_TEST_COUNT
    6. Run: python check_raptor.py

Sections:
     1. Project Structure        9.  Test Suite
     2. setup.py Discovery      10.  CLI Entry Point
     3. Core Module Imports      11.  Requirements & Environments
     4. Acquisition Imports      12.  Code Quality
     5. Class Methods            13.  Git Repository
     6. Dependencies             14.  Cache & Data Directories
     7. Version Consistency      15.  .gitignore
     8. Dashboard                16.  Smoke Tests

Author: Ayeh Bolouki
Version: 2.2.2
"""

import sys
import os
import re
import importlib
import subprocess
import time
import traceback
from pathlib import Path
from io import StringIO

# ─────────────────────────────────────────────────────────────────────
# Configuration — UPDATE WHEN ADDING NEW MODULES
# ─────────────────────────────────────────────────────────────────────

# -- Packages that must be discoverable by find_packages() --
EXPECTED_PACKAGES = [
    'raptor',
    'raptor.utils',
    'raptor.pipelines',
    'raptor.external_modules',
    'raptor.external_modules.module6_de_analysis',
    'raptor.external_modules.acquisition',
    'raptor.dashboard',
    'raptor.dashboard.components',
    'raptor.dashboard.pages',
]

# -- Core RAPTOR module files (non-optional, must import) --
EXPECTED_CORE_MODULES = [
    'raptor.simulation',
    'raptor.profiler',
    'raptor.quality_assessment',
    'raptor.recommender',
    'raptor.ml_recommender',
    'raptor.synthetic_benchmarks',
]

# -- Optional RAPTOR module files (may fail if deps missing) --
EXPECTED_OPTIONAL_MODULES = [
    ('raptor.de_import',              'Module 7: DE Import'),
    ('raptor.parameter_optimization', 'Module 8: Parameter Optimization'),
    ('raptor.ensemble',               'Module 9: Ensemble Analysis'),
]

# -- Acquisition subpackage modules --
EXPECTED_ACQUISITION_MODULES = [
    'raptor.external_modules.acquisition.base',
    'raptor.external_modules.acquisition.datasets',
    'raptor.external_modules.acquisition.cache',
    'raptor.external_modules.acquisition.catalog',
    'raptor.external_modules.acquisition.geo',
    'raptor.external_modules.acquisition.tcga',
    'raptor.external_modules.acquisition.arrayexpress',
    'raptor.external_modules.acquisition.sra',
    'raptor.external_modules.acquisition.gene_mapping',
    'raptor.external_modules.acquisition.pooling',
]

# -- Utils module files --
EXPECTED_UTILS_MODULES = [
    'raptor.utils.validation',
    'raptor.utils.errors',
    'raptor.utils.sample_sheet',
]

# -- Key exports from raptor/__init__.py (critical subset) --
EXPECTED_RAPTOR_EXPORTS = [
    # Core (always available)
    'SimulationConfig', 'RNAseqSimulator', 'simulate_rnaseq',
    'DataProfile', 'RNAseqDataProfiler', 'profile_data_quick',
    'DataQualityAssessor', 'quick_quality_check',
    'PipelineRecommender', 'recommend_pipeline',
    'MLPipelineRecommender', 'train_ml_recommender',
    'SyntheticBenchmarkGenerator', 'generate_training_data',
    # Utils
    'RAPTORError', 'ValidationError',
    # Optional (may be None but must exist as attributes)
    'DEResult', 'import_de_results',
    'ParameterOptimizer', 'optimize_with_ground_truth', 'optimize_with_fdr_control',
    'EnsembleResult', 'ensemble_fisher', 'ensemble_brown',
    # Acquisition (may be None)
    'AcquiredDataset', 'GEOConnector', 'SRAConnector', 'GeneIDMapper',
    # Helper functions
    'get_version', 'get_info', 'get_available_modules', 'validate_installation',
]

# -- Acquisition subpackage exports --
EXPECTED_ACQUISITION_EXPORTS = [
    'AcquiredDataset', 'PooledDataset', 'CacheManager', 'DataCatalog',
    'GEOConnector', 'TCGAConnector', 'ArrayExpConnector', 'SRAConnector',
    'GeneIDMapper', 'PoolingEngine', 'SearchResult',
    'get_available_components', 'DEFAULT_CACHE_DIR',
]

# -- Key methods that must exist on specific classes --
EXPECTED_CLASS_METHODS = {
    'SRAConnector': [
        '_extract_gsm_ids', '_gse_from_gsm', 'find_linked_gse',
        'search', 'download', 'get_run_table', 'generate_download_commands',
    ],
    'GEOConnector': [
        'search', 'download', 'get_study_info', '_detect_gene_id_type',
    ],
    'TCGAConnector': [
        'search', 'download',
    ],
    'GeneIDMapper': [
        'convert', 'convert_index', 'detect_id_type', 'harmonize_to_common',
    ],
    'PoolingEngine': [
        'merge',
    ],
    'DataQualityAssessor': [],   # discovered at runtime
    'RNAseqDataProfiler': [],    # discovered at runtime
}

# -- Dependencies --
REQUIRED_DEPS = [
    'numpy', 'pandas', 'scipy', 'click', 'sklearn', 'statsmodels',
    'matplotlib', 'seaborn', 'plotly', 'yaml', 'tqdm', 'joblib',
    'toml', 'colorama',
]

OPTIONAL_DEPS = {
    'pyarrow':    ('Parquet caching (acquisition)', 'pyarrow'),
    'requests':   ('GEO/SRA/TCGA/ArrayExpress API access', 'requests'),
    'GEOparse':   ('GEO dataset parsing', 'GEOparse'),
    'Bio':        ('SRA/GEO Entrez search', 'biopython'),
    'mygene':     ('Gene ID mapping (MyGene.info)', 'mygene'),
    'streamlit':  ('Dashboard', 'streamlit'),
    'combat':     ('ComBat batch correction (pooling)', 'combat'),
}

# -- Dashboard --
EXPECTED_DASHBOARD_PAGES = [
    '01_', '02_', '03_', '04_', '05_',
    '06_', '07_', '08_', '09_', '10_',
]

# -- Tests --
EXPECTED_TEST_COUNT = 105
EXPECTED_TEST_DIR = 'tests'

# -- Version files --
VERSION_LOCATIONS = [
    ('raptor/__init__.py', r'__version__\s*=\s*["\x27]([^"\x27]+)["\x27]'),
    ('setup.py',           r'VERSION\s*=\s*["\x27]([^"\x27]+)["\x27]'),
]

# -- Code directories to scan --
CODE_DIRECTORIES = [
    'raptor/',
    'tests/',
]

# -- Hardcoded version patterns to flag --
HARDCODED_VERSION_PATTERNS = [
    (r'Version:\s*2\.\d+\.\d+', 'docstring version'),
    (r"version='2\.\d+\.\d+'", 'hardcoded version string'),
]

# -- Core files that must exist --
EXPECTED_ROOT_FILES = [
    'setup.py', 'requirements.txt', 'README.md', 'CHANGELOG.md',
    '.gitignore',
]

EXPECTED_CORE_FILES = [
    'raptor/__init__.py',
    'raptor/cli.py',
    'raptor/simulation.py',
    'raptor/profiler.py',
    'raptor/quality_assessment.py',
    'raptor/recommender.py',
    'raptor/ml_recommender.py',
    'raptor/synthetic_benchmarks.py',
    'raptor/de_import.py',
    'raptor/parameter_optimization.py',
    'raptor/ensemble.py',
    'raptor/utils/__init__.py',
    'raptor/utils/validation.py',
    'raptor/utils/errors.py',
    'raptor/utils/sample_sheet.py',
    'raptor/external_modules/__init__.py',
    'raptor/external_modules/acquisition/__init__.py',
    'raptor/dashboard/app.py',
    'raptor/dashboard/components/sidebar.py',
]


# ─────────────────────────────────────────────────────────────────────
# Display Helpers
# ─────────────────────────────────────────────────────────────────────

if sys.platform == 'win32':
    os.system('')  # Enable ANSI on Windows

PASS = "\033[92m PASS \033[0m"
FAIL = "\033[91m FAIL \033[0m"
WARN = "\033[93m WARN \033[0m"
INFO = "\033[94m INFO \033[0m"
SKIP = "\033[90m SKIP \033[0m"

results = {'pass': 0, 'fail': 0, 'warn': 0, 'skip': 0}
all_output = StringIO()
verbose = '--verbose' in sys.argv
show_fix = '--fix' in sys.argv


def log(msg):
    print(msg)
    all_output.write(msg + '\n')


def check(condition, message, fix_hint="", warn_only=False):
    if condition:
        log(f"  [{PASS}] {message}")
        results['pass'] += 1
    elif warn_only:
        log(f"  [{WARN}] {message}")
        if fix_hint and show_fix:
            log(f"          Fix: {fix_hint}")
        results['warn'] += 1
    else:
        log(f"  [{FAIL}] {message}")
        if fix_hint and show_fix:
            log(f"          Fix: {fix_hint}")
        results['fail'] += 1


def section(number, title):
    log(f"\n{'='*65}")
    log(f"  {number}. {title}")
    log(f"{'='*65}")


def should_run(section_num):
    for arg in sys.argv:
        if arg.startswith('--section'):
            idx = sys.argv.index(arg)
            if idx + 1 < len(sys.argv):
                sections = sys.argv[idx + 1].split(',')
                return str(section_num) in sections
            return True
    return True


def is_quick():
    return '--quick' in sys.argv


# ─────────────────────────────────────────────────────────────────────
# 1. Project Structure
# ─────────────────────────────────────────────────────────────────────

def check_project_structure():
    section(1, "Project Structure")
    root = Path('.')

    for f in EXPECTED_ROOT_FILES:
        check((root / f).exists(), f"{f} exists",
              f"Missing root file: {f}")

    for f in EXPECTED_CORE_FILES:
        check((root / f).exists(), f"{f}",
              f"Missing: {f}")

    acq_dir = root / 'raptor' / 'external_modules' / 'acquisition'
    expected_acq = [
        'base.py', 'datasets.py', 'cache.py', 'catalog.py',
        'geo.py', 'tcga.py', 'arrayexpress.py', 'sra.py',
        'gene_mapping.py', 'pooling.py',
    ]
    for f in expected_acq:
        check((acq_dir / f).exists(), f"acquisition/{f}",
              f"Missing: raptor/external_modules/acquisition/{f}")

    tests_dir = root / 'tests'
    if tests_dir.exists():
        test_files = list(tests_dir.glob('test_*.py'))
        log(f"  [{INFO}] Found {len(test_files)} test files")
        check((tests_dir / 'test_acquisition.py').exists(), "test_acquisition.py exists")
    else:
        check(False, "tests/ directory exists", "mkdir tests")

    r_scripts = root / 'raptor' / 'external_modules' / 'module6_de_analysis' / 'r_scripts'
    if r_scripts.exists():
        r_files = list(r_scripts.glob('*.R'))
        log(f"  [{INFO}] {len(r_files)} R scripts in module6_de_analysis")


# ─────────────────────────────────────────────────────────────────────
# 2. setup.py Package Discovery
# ─────────────────────────────────────────────────────────────────────

def check_setup_py():
    section(2, "setup.py Package Discovery")
    setup_text = Path('setup.py').read_text(encoding='utf-8')
    uses_find_packages = 'find_packages' in setup_text

    if uses_find_packages:
        log(f"  [{INFO}] setup.py uses find_packages() — auto-discovers subpackages")
        for pkg in EXPECTED_PACKAGES:
            pkg_path = Path(pkg.replace('.', '/'))
            check((pkg_path / '__init__.py').exists(),
                  f"Package '{pkg}' discoverable ({pkg_path}/__init__.py)",
                  f"Create {pkg_path}/__init__.py so find_packages() discovers it")
    else:
        log(f"  [{INFO}] setup.py uses explicit package list")
        for pkg in EXPECTED_PACKAGES:
            check(pkg in setup_text, f"setup.py lists '{pkg}'",
                  f"Add '{pkg}' to packages list in setup.py")

    if 'entry_points' in setup_text or 'console_scripts' in setup_text:
        check(True, "CLI entry point configured")
    else:
        check(False, "CLI entry point (console_scripts) configured",
              "Add entry_points to setup.py", warn_only=True)

    for extra in ['dashboard', 'acquisition', 'dev', 'all']:
        check(f'"{extra}"' in setup_text or f"'{extra}'" in setup_text,
              f"extras_require['{extra}'] defined",
              f"Add '{extra}' extras to setup.py", warn_only=True)

    for dep in ['requests', 'pyarrow', 'click', 'numpy', 'pandas']:
        check(dep in setup_text, f"setup.py includes {dep}")


# ─────────────────────────────────────────────────────────────────────
# 3. Core Module Imports
# ─────────────────────────────────────────────────────────────────────

def check_imports():
    section(3, "Core Module Imports")

    try:
        import raptor
        importlib.reload(raptor)
        check(True, f"import raptor (v{raptor.__version__})")
    except Exception as e:
        check(False, f"import raptor — {e}")
        return

    for mod_path in EXPECTED_CORE_MODULES:
        mod_name = mod_path.split('.')[-1]
        try:
            importlib.import_module(mod_path)
            check(True, f"import {mod_name}")
        except Exception as e:
            check(False, f"import {mod_name} — {e}")

    for mod_path, description in EXPECTED_OPTIONAL_MODULES:
        mod_name = mod_path.split('.')[-1]
        try:
            importlib.import_module(mod_path)
            check(True, f"import {mod_name} ({description})")
        except Exception as e:
            check(False, f"import {mod_name} — {e}", warn_only=True)

    for mod_path in EXPECTED_UTILS_MODULES:
        mod_name = mod_path.split('.')[-1]
        try:
            importlib.import_module(mod_path)
            check(True, f"import utils.{mod_name}")
        except Exception as e:
            check(False, f"import utils.{mod_name} — {e}")

    # raptor.__init__ exports
    missing_exports = []
    for name in EXPECTED_RAPTOR_EXPORTS:
        if not hasattr(raptor, name):
            missing_exports.append(name)

    check(len(missing_exports) == 0,
          f"All {len(EXPECTED_RAPTOR_EXPORTS)} key raptor exports available",
          f"Missing from raptor/__init__.py: {missing_exports}")

    # get_available_modules()
    try:
        modules = raptor.get_available_modules()
        available = [k for k, v in modules.items() if v]
        unavailable = [k for k, v in modules.items() if not v]
        log(f"  [{INFO}] get_available_modules(): {len(available)} available, "
            f"{len(unavailable)} unavailable")
        if unavailable:
            log(f"  [{INFO}] Unavailable: {', '.join(unavailable)}")
    except Exception as e:
        check(False, f"get_available_modules() — {e}", warn_only=True)

    # validate_installation()
    try:
        report = raptor.validate_installation()
        check(report['exports_valid'],
              "validate_installation(): all exports valid",
              f"Missing exports: {report.get('missing_exports', [])[:5]}")
    except Exception as e:
        check(False, f"validate_installation() — {e}", warn_only=True)


# ─────────────────────────────────────────────────────────────────────
# 4. Acquisition Module Imports
# ─────────────────────────────────────────────────────────────────────

def check_acquisition_imports():
    section(4, "Acquisition Module Imports")

    try:
        from raptor.external_modules import acquisition
        check(True, "from raptor.external_modules import acquisition")
    except Exception as e:
        check(False, f"acquisition import — {e}",
              "Check raptor/external_modules/__init__.py")
        return

    for mod_path in EXPECTED_ACQUISITION_MODULES:
        mod_name = mod_path.split('.')[-1]
        try:
            importlib.import_module(mod_path)
            check(True, f"import {mod_name}")
        except Exception as e:
            check(False, f"import {mod_name} — {e}")

    missing_exports = []
    for name in EXPECTED_ACQUISITION_EXPORTS:
        try:
            obj = getattr(acquisition, name)
            if obj is None:
                missing_exports.append(f"{name} (is None)")
        except AttributeError:
            missing_exports.append(name)

    check(len(missing_exports) == 0,
          f"All {len(EXPECTED_ACQUISITION_EXPORTS)} acquisition exports available",
          f"Missing: {missing_exports}" if missing_exports else "")

    try:
        for mod_path in EXPECTED_ACQUISITION_MODULES:
            importlib.reload(importlib.import_module(mod_path))
        check(True, "No circular imports (all acquisition modules reloadable)")
    except Exception as e:
        check(False, f"Circular import detected: {e}", warn_only=True)


# ─────────────────────────────────────────────────────────────────────
# 5. Class Methods Verification
# ─────────────────────────────────────────────────────────────────────

def check_class_methods():
    section(5, "Class Methods & Signatures")

    class_map = {}

    try:
        from raptor.external_modules.acquisition import (
            SRAConnector, GEOConnector, TCGAConnector,
            GeneIDMapper, PoolingEngine,
        )
        class_map.update({
            'SRAConnector': SRAConnector,
            'GEOConnector': GEOConnector,
            'TCGAConnector': TCGAConnector,
            'GeneIDMapper': GeneIDMapper,
            'PoolingEngine': PoolingEngine,
        })
    except ImportError as e:
        log(f"  [{WARN}] Cannot import acquisition classes: {e}")

    try:
        from raptor import DataQualityAssessor, RNAseqDataProfiler
        class_map.update({
            'DataQualityAssessor': DataQualityAssessor,
            'RNAseqDataProfiler': RNAseqDataProfiler,
        })
    except ImportError as e:
        log(f"  [{WARN}] Cannot import core classes: {e}")

    for cls_name, expected_methods in EXPECTED_CLASS_METHODS.items():
        cls = class_map.get(cls_name)
        if cls is None:
            check(False, f"{cls_name} class not available")
            continue

        if not expected_methods:
            # Discovery mode: list public methods for classes we haven't verified yet
            public_methods = [m for m in dir(cls)
                              if not m.startswith('_') and callable(getattr(cls, m, None))]
            log(f"  [{INFO}] {cls_name} has {len(public_methods)} public methods: "
                f"{', '.join(public_methods[:8])}{'...' if len(public_methods) > 8 else ''}")
            check(len(public_methods) > 0, f"{cls_name} has public methods")
            continue

        for method_name in expected_methods:
            has_method = hasattr(cls, method_name)
            check(has_method, f"{cls_name}.{method_name}()",
                  f"Add {method_name} method to {cls_name}")

            if has_method and verbose:
                import inspect
                try:
                    sig = inspect.signature(getattr(cls, method_name))
                    log(f"          Signature: {method_name}{sig}")
                except (ValueError, TypeError):
                    pass


# ─────────────────────────────────────────────────────────────────────
# 6. Dependencies
# ─────────────────────────────────────────────────────────────────────

def check_dependencies():
    section(6, "Dependencies")

    installed = 0
    missing_count = 0
    for dep in REQUIRED_DEPS:
        try:
            mod = importlib.import_module(dep)
            version = getattr(mod, '__version__', '?')
            check(True, f"{dep} {version}")
            installed += 1
        except ImportError:
            check(False, f"{dep} — NOT INSTALLED", f"pip install {dep}")
            missing_count += 1

    log(f"  [{INFO}] Required: {installed} installed, {missing_count} missing")

    for dep, (purpose, pip_name) in OPTIONAL_DEPS.items():
        try:
            mod = importlib.import_module(dep)
            version = getattr(mod, '__version__', '?')
            check(True, f"{dep} {version} — {purpose}")
        except ImportError:
            check(False, f"{dep} — {purpose} — NOT INSTALLED",
                  f"pip install {pip_name}", warn_only=True)


# ─────────────────────────────────────────────────────────────────────
# 7. Version Consistency
# ─────────────────────────────────────────────────────────────────────

def check_versions():
    section(7, "Version Consistency")
    versions = {}

    for filepath, pattern in VERSION_LOCATIONS:
        try:
            text = Path(filepath).read_text(encoding='utf-8')
            match = re.search(pattern, text)
            if match:
                versions[filepath] = match.group(1)
                log(f"  [{INFO}] {filepath}: {match.group(1)}")
            elif filepath == 'setup.py' and '__version__' in text:
                log(f"  [{INFO}] {filepath}: reads version dynamically (good practice)")
            else:
                log(f"  [{WARN}] {filepath}: version pattern not found")
        except FileNotFoundError:
            log(f"  [{FAIL}] {filepath}: file not found")

    try:
        from raptor import cli
        if hasattr(cli, '__version__'):
            versions['cli.py'] = cli.__version__
            log(f"  [{INFO}] cli.py: {cli.__version__}")
    except Exception:
        pass

    try:
        from raptor.external_modules import acquisition
        if hasattr(acquisition, '__version__'):
            acq_ver = acquisition.__version__
            log(f"  [{INFO}] acquisition/__init__.py: {acq_ver}")
            import raptor
            if acq_ver != raptor.__version__:
                check(False,
                      f"Acquisition version ({acq_ver}) != raptor ({raptor.__version__})",
                      "Update __version__ in acquisition/__init__.py",
                      warn_only=True)
    except Exception:
        pass

    unique_versions = set(versions.values())
    if len(unique_versions) == 1:
        check(True, f"All versions match: {unique_versions.pop()}")
    elif len(unique_versions) > 1:
        check(False, f"Version MISMATCH: {versions}",
              "Update all files to the same version")
    else:
        check(False, "Could not determine any version")

    hardcoded_count = 0
    for dirpath in CODE_DIRECTORIES:
        for pyfile in Path(dirpath).rglob('*.py'):
            try:
                content = pyfile.read_text(encoding='utf-8')
                for pattern, desc in HARDCODED_VERSION_PATTERNS:
                    hardcoded_count += len(re.findall(pattern, content))
            except Exception:
                pass
    if hardcoded_count > 0:
        log(f"  [{INFO}] {hardcoded_count} hardcoded version references in codebase")


# ─────────────────────────────────────────────────────────────────────
# 8. Dashboard
# ─────────────────────────────────────────────────────────────────────

def check_dashboard():
    section(8, "Dashboard")
    pages_dir = Path('raptor/dashboard/pages')
    if not pages_dir.exists():
        check(False, "Dashboard pages directory exists")
        return

    page_files = sorted([
        f.name for f in pages_dir.glob('*.py') if not f.name.startswith('__')
    ])
    log(f"  [{INFO}] Found {len(page_files)} page files")

    for prefix in EXPECTED_DASHBOARD_PAGES:
        matching = [f for f in page_files if f.startswith(prefix)]
        check(len(matching) == 1,
              f"Page {prefix}* ({matching[0] if matching else 'MISSING'})",
              f"Missing dashboard page with prefix {prefix}")

    prefixes = [f[:3] for f in page_files]
    duplicates = [p for p in set(prefixes) if prefixes.count(p) > 1]
    check(len(duplicates) == 0, "No duplicate page numbers",
          f"Duplicate prefixes: {duplicates}")

    syntax_errors = []
    for f in page_files:
        try:
            import py_compile
            py_compile.compile(str(pages_dir / f), doraise=True)
        except py_compile.PyCompileError as e:
            syntax_errors.append(f"{f}: {e}")

    check(len(syntax_errors) == 0,
          f"All {len(page_files)} pages have valid syntax",
          f"Syntax errors: {syntax_errors}")

    app_file = Path('raptor/dashboard/app.py')
    if app_file.exists():
        app_text = app_file.read_text(encoding='utf-8')
        check('set_page_config' in app_text, "app.py has st.set_page_config()")
        check('Data Acquisition' in app_text or 'acquisition' in app_text.lower(),
              "app.py references Data Acquisition",
              "Update app.py sidebar to include Data Acquisition", warn_only=True)

    sidebar_file = Path('raptor/dashboard/components/sidebar.py')
    check(sidebar_file.exists(), "components/sidebar.py exists")


# ─────────────────────────────────────────────────────────────────────
# 9. Test Suite
# ─────────────────────────────────────────────────────────────────────

def check_tests():
    section(9, "Test Suite")
    if is_quick():
        log(f"  [{SKIP}] Skipped (--quick mode)")
        results['skip'] += 1
        return

    tests_dir = Path(EXPECTED_TEST_DIR)
    if not tests_dir.exists():
        check(False, "tests/ directory exists")
        return

    test_files = sorted(tests_dir.glob('test_*.py'))
    log(f"  [{INFO}] Test files: {', '.join(f.name for f in test_files)}")

    try:
        result = subprocess.run(
            [sys.executable, '-m', 'pytest', str(tests_dir),
             '--collect-only', '-q'],
            capture_output=True, text=True, timeout=30,
        )
        match = re.search(r'(\d+) test', result.stdout)
        if match:
            n = int(match.group(1))
            check(n >= EXPECTED_TEST_COUNT,
                  f"Collected {n} tests (expected >= {EXPECTED_TEST_COUNT})")
        else:
            check(False, "Could not parse test count")
            if verbose:
                log(f"          stdout: {result.stdout[:200]}")
    except Exception as e:
        check(False, f"Test collection failed: {e}")

    log(f"  [{INFO}] Running full test suite...")
    try:
        result = subprocess.run(
            [sys.executable, '-m', 'pytest', str(tests_dir),
             '-v', '--tb=short', '-q'],
            capture_output=True, text=True, timeout=600,
        )
        passed = int(m.group(1)) if (m := re.search(r'(\d+) passed', result.stdout)) else 0
        failed = int(m.group(1)) if (m := re.search(r'(\d+) failed', result.stdout)) else 0
        errors = int(m.group(1)) if (m := re.search(r'(\d+) error', result.stdout)) else 0

        check(failed == 0 and errors == 0,
              f"All tests pass ({passed} passed, {failed} failed, {errors} errors)",
              "Run: pytest tests/ -v --tb=long")

        if (failed > 0 or errors > 0) and verbose:
            for line in result.stdout.split('\n'):
                if 'FAILED' in line or 'ERROR' in line:
                    log(f"          {line.strip()}")

        time_match = re.search(r'in ([\d.]+)s', result.stdout)
        if time_match:
            log(f"  [{INFO}] Completed in {time_match.group(1)}s")

    except subprocess.TimeoutExpired:
        check(False, "Test run timed out (>600s)")
    except Exception as e:
        check(False, f"Test run failed: {e}")


# ─────────────────────────────────────────────────────────────────────
# 10. CLI Entry Point
# ─────────────────────────────────────────────────────────────────────

def check_cli():
    section(10, "CLI Entry Point")
    if is_quick():
        log(f"  [{SKIP}] Skipped (--quick mode)")
        results['skip'] += 1
        return

    try:
        result = subprocess.run(
            [sys.executable, '-m', 'raptor.cli', '--version'],
            capture_output=True, text=True, timeout=15,
        )
        if result.returncode == 0:
            cli_output = result.stdout.strip()
            check(True, f"raptor --version: {cli_output}")

            try:
                import raptor
                init_ver = raptor.__version__
                if init_ver not in cli_output:
                    check(False,
                          f"CLI version mismatch: CLI='{cli_output}', __init__='{init_ver}'",
                          "Use: @click.version_option(version=raptor.__version__)",
                          warn_only=True)
            except Exception:
                pass
        else:
            check(False, "raptor --version failed",
                  "Check CLI entry point in setup.py", warn_only=True)
            if verbose:
                log(f"          stderr: {result.stderr[:200]}")
    except Exception as e:
        check(False, f"CLI check failed: {e}", warn_only=True)

    try:
        result = subprocess.run(
            [sys.executable, '-m', 'raptor.cli', '--help'],
            capture_output=True, text=True, timeout=15,
        )
        if result.returncode == 0:
            commands = [c for c in ['quantify', 'dashboard', 'de', 'optimize', 'ensemble']
                        if c in result.stdout]
            check(len(commands) >= 3,
                  f"CLI help lists {len(commands)} commands: {', '.join(commands)}")
    except Exception as e:
        check(False, f"CLI help failed: {e}", warn_only=True)


# ─────────────────────────────────────────────────────────────────────
# 11. Requirements & Environment Files
# ─────────────────────────────────────────────────────────────────────

def check_requirements():
    section(11, "Requirements & Environment Files")
    req_file = Path('requirements.txt')
    if not req_file.exists():
        check(False, "requirements.txt exists")
        return

    req_text = req_file.read_text(encoding='utf-8').lower()
    for dep in ['numpy', 'pandas', 'click', 'scipy', 'scikit-learn',
                'statsmodels', 'matplotlib', 'plotly', 'pyyaml']:
        check(dep in req_text, f"requirements.txt includes {dep}")

    for dep in ['requests', 'pyarrow']:
        check(dep in req_text, f"requirements.txt includes {dep} (acquisition)")

    for env_file in ['environment.yml', 'environment-full.yml']:
        path = Path(env_file)
        if path.exists():
            env_text = path.read_text(encoding='utf-8').lower()
            has_acq = 'acquisition' in env_text or 'requests' in env_text
            log(f"  [{INFO}] {env_file} ({path.stat().st_size:,} bytes)"
                f"{' — has acquisition deps' if has_acq else ''}")


# ─────────────────────────────────────────────────────────────────────
# 12. Code Quality
# ─────────────────────────────────────────────────────────────────────

def check_code_quality():
    section(12, "Code Quality")

    encoding_issues = []
    total_files = 0
    for dirpath in CODE_DIRECTORIES:
        for pyfile in Path(dirpath).rglob('*.py'):
            total_files += 1
            try:
                pyfile.read_text(encoding='utf-8')
            except UnicodeDecodeError:
                encoding_issues.append(str(pyfile))

    check(len(encoding_issues) == 0, f"All {total_files} Python files are valid UTF-8",
          f"Encoding issues in: {encoding_issues}")

    # print() in library code (skip dashboard pages)
    print_locs = []
    for pyfile in Path('raptor').rglob('*.py'):
        if 'dashboard' in str(pyfile):
            continue
        try:
            lines = pyfile.read_text(encoding='utf-8').split('\n')
            in_main_guard = False
            for i, line in enumerate(lines, 1):
                stripped = line.strip()
                if "if __name__" in stripped:
                    in_main_guard = True
                if in_main_guard:
                    continue
                if stripped.startswith('print(') and not stripped.startswith('#'):
                    print_locs.append(f"{pyfile.relative_to('.')}:{i}")
        except Exception:
            pass

    if print_locs:
        check(False, f"print() in library code ({len(print_locs)} locations)",
              f"Use logging instead. Found in: {print_locs[:5]}", warn_only=True)
    else:
        check(True, "No print() in library code (uses logging)")

    syntax_errors = []
    for dirpath in CODE_DIRECTORIES:
        for pyfile in Path(dirpath).rglob('*.py'):
            try:
                import py_compile
                py_compile.compile(str(pyfile), doraise=True)
            except py_compile.PyCompileError as e:
                syntax_errors.append(str(e)[:100])

    check(len(syntax_errors) == 0, f"All {total_files} Python files have valid syntax",
          f"Syntax errors: {syntax_errors}")

    missing_enc = []
    for pyfile in Path('raptor').rglob('*.py'):
        try:
            content = pyfile.read_text(encoding='utf-8')
            for i, line in enumerate(content.split('\n'), 1):
                if ('open(' in line and ("'w'" in line or '"w"' in line)
                        and 'encoding' not in line and not line.strip().startswith('#')):
                    missing_enc.append(f"{pyfile.relative_to('.')}:{i}")
        except Exception:
            pass

    if missing_enc:
        check(False, f"Missing encoding='utf-8' in {len(missing_enc)} file writes",
              f"Found in: {missing_enc[:5]}", warn_only=True)
    else:
        check(True, "All file writes use encoding='utf-8'")


# ─────────────────────────────────────────────────────────────────────
# 13. Git Repository
# ─────────────────────────────────────────────────────────────────────

def check_git():
    section(13, "Git Repository")
    try:
        result = subprocess.run(['git', 'status', '--porcelain'],
                                capture_output=True, text=True, timeout=15)
        if result.returncode != 0:
            check(False, "Git repository accessible")
            return

        lines = [l for l in result.stdout.strip().split('\n') if l.strip()]
        modified = [l for l in lines if l.startswith(' M') or l.startswith('M ')]
        untracked = [l for l in lines if l.startswith('??')]

        if not lines:
            check(True, "Working directory clean")
        else:
            check(len(modified) == 0,
                  f"Uncommitted: {len(modified)} modified, {len(untracked)} untracked",
                  "git add and commit your changes", warn_only=True)
            if verbose and modified:
                for l in modified[:10]:
                    log(f"          {l.strip()}")

        branch = subprocess.run(['git', 'branch', '--show-current'],
                                capture_output=True, text=True, timeout=10)
        if branch.returncode == 0:
            log(f"  [{INFO}] Branch: {branch.stdout.strip()}")

        status = subprocess.run(['git', 'status', '-sb'],
                                capture_output=True, text=True, timeout=10)
        if 'ahead' in status.stdout:
            check(False, "Local ahead of remote — push needed",
                  "git push origin main", warn_only=True)
        elif 'behind' in status.stdout:
            check(False, "Local behind remote — pull needed",
                  "git pull origin main", warn_only=True)
        else:
            check(True, "In sync with remote")

    except FileNotFoundError:
        check(False, "git not on PATH", warn_only=True)
    except Exception as e:
        check(False, f"Git check failed: {e}", warn_only=True)


# ─────────────────────────────────────────────────────────────────────
# 14. Cache & Data Directories
# ─────────────────────────────────────────────────────────────────────

def check_cache():
    section(14, "Cache & Data Directories")
    cache_dir = Path.home() / '.raptor_cache'
    if cache_dir.exists():
        n_files = sum(1 for _ in cache_dir.rglob('*') if _.is_file())
        total_size = sum(f.stat().st_size for f in cache_dir.rglob('*') if f.is_file())
        log(f"  [{INFO}] Cache: {cache_dir} ({n_files} files, {total_size / 1024 / 1024:.1f} MB)")
    else:
        log(f"  [{INFO}] No cache directory yet (created on first download)")

    parquet_dir = Path.home() / '.raptor' / 'parquet_cache'
    if parquet_dir.exists():
        n_pq = sum(1 for _ in parquet_dir.rglob('*.parquet'))
        log(f"  [{INFO}] Parquet cache: {n_pq} files in {parquet_dir}")

    pycache_dirs = list(Path('raptor').rglob('__pycache__'))
    log(f"  [{INFO}] {len(pycache_dirs)} __pycache__ directories")
    check(True, "Cache check complete")


# ─────────────────────────────────────────────────────────────────────
# 15. .gitignore
# ─────────────────────────────────────────────────────────────────────

def check_gitignore():
    section(15, ".gitignore")
    gi_file = Path('.gitignore')
    if not gi_file.exists():
        check(False, ".gitignore exists")
        return

    gi_text = gi_file.read_text(encoding='utf-8')

    for pattern, desc in [
        ('.venv', 'Virtual environments'),
        ('__pycache__', 'Bytecode cache'),
        ('*.egg-info', 'Package metadata'),
        ('.pytest_cache', 'Pytest cache'),
    ]:
        check(pattern in gi_text, f".gitignore includes {pattern} ({desc})")

    for pattern, desc in [
        ('.raptor_cache', 'RAPTOR data cache'),
        ('dist/', 'Build distributions'),
        ('build/', 'Build directory'),
    ]:
        if pattern not in gi_text:
            check(False, f".gitignore includes {pattern} ({desc})",
                  f"Add {pattern} to .gitignore", warn_only=True)
        else:
            check(True, f".gitignore includes {pattern} ({desc})")


# ─────────────────────────────────────────────────────────────────────
# 16. Smoke Tests (offline functionality)
# ─────────────────────────────────────────────────────────────────────

def check_smoke_test():
    section(16, "Smoke Tests (offline functionality)")
    if is_quick():
        log(f"  [{SKIP}] Skipped (--quick mode)")
        results['skip'] += 1
        return

    try:
        import numpy as np
        import pandas as pd

        # ── Core module smoke tests ──
        log(f"  [{INFO}] Testing core modules...")

        np.random.seed(42)
        counts = pd.DataFrame(
            np.random.negative_binomial(5, 0.3, (500, 6)),
            index=[f'GENE_{i}' for i in range(500)],
            columns=['ctrl_1', 'ctrl_2', 'ctrl_3', 'treat_1', 'treat_2', 'treat_3'],
        )
        counts.index.name = 'gene_id'
        meta = pd.DataFrame(
            {'condition': ['ctrl', 'ctrl', 'ctrl', 'treat', 'treat', 'treat']},
            index=counts.columns,
        )
        meta.index.name = 'sample_id'

        # Simulation
        try:
            from raptor import simulate_rnaseq
            sim_result = simulate_rnaseq()
            check(sim_result is not None and hasattr(sim_result, 'counts'),
                  "simulate_rnaseq() produces output")
        except Exception as e:
            check(False, f"simulate_rnaseq() — {e}")

        # Profiler
        try:
            from raptor import profile_data_quick
            # profile_data_quick may expect sample_id as column, not index
            meta_col = meta.copy()
            if 'sample_id' not in meta_col.columns:
                meta_col = meta_col.reset_index()
            profile = profile_data_quick(counts, meta_col)
            check(profile is not None,
                  f"profile_data_quick() — BCV={getattr(profile, 'bcv', '?')}")
        except Exception as e:
            check(False, f"profile_data_quick() — {e}")

        # Quality Assessment
        try:
            from raptor import DataQualityAssessor
            import inspect
            sig = inspect.signature(DataQualityAssessor.__init__)
            params = list(sig.parameters.keys())
            log(f"  [{INFO}] DataQualityAssessor.__init__ params: {params[1:]}")
            # Try constructing with detected signature
            meta_col = meta.copy()
            if 'sample_id' not in meta_col.columns:
                meta_col = meta_col.reset_index()
            if 'group_column' in params:
                assessor = DataQualityAssessor(counts, meta_col, group_column='condition')
            elif 'metadata' in params:
                assessor = DataQualityAssessor(counts, meta_col)
            else:
                assessor = DataQualityAssessor(counts)
            check(assessor is not None, "DataQualityAssessor created successfully")
        except Exception as e:
            check(False, f"DataQualityAssessor — {e}")

        # Recommender
        try:
            from raptor import recommend_pipeline
            import inspect
            sig = inspect.signature(recommend_pipeline)
            n_params = len([p for p in sig.parameters.values()
                           if p.default is inspect.Parameter.empty])
            if n_params >= 2:
                rec = recommend_pipeline(counts, meta)
            else:
                rec = recommend_pipeline(counts)
            check(rec is not None,
                  f"recommend_pipeline() — primary: {getattr(rec, 'primary_pipeline', '?')}")
        except Exception as e:
            check(False, f"recommend_pipeline() — {e}")

        # ── Acquisition smoke tests ──
        log(f"  [{INFO}] Testing acquisition modules...")

        from raptor.external_modules.acquisition import (
            AcquiredDataset, PooledDataset, GeneIDMapper,
            PoolingEngine, SRAConnector,
        )

        acq_counts = pd.DataFrame(
            np.random.randint(0, 100, (200, 4)),
            index=[f'GENE_{i}' for i in range(200)],
            columns=['S1', 'S2', 'S3', 'S4'],
        )
        acq_counts.index.name = 'gene_id'
        acq_meta = pd.DataFrame(
            {'condition': ['ctrl', 'ctrl', 'treat', 'treat']},
            index=['S1', 'S2', 'S3', 'S4'],
        )
        acq_meta.index.name = 'sample_id'

        ds = AcquiredDataset(
            counts_df=acq_counts, metadata=acq_meta,
            source_info={'repository': 'test', 'accession': 'TEST001'},
            gene_id_type='symbol',
        )
        check(ds.n_genes == 200 and ds.n_samples == 4, "AcquiredDataset creation")

        integrity = ds.validate_integrity()
        check(integrity['valid'], "AcquiredDataset.validate_integrity()")

        mapper = GeneIDMapper('Homo sapiens')
        check(mapper.detect_id_type(['TP53', 'BRCA1', 'EGFR']) == 'symbol',
              "GeneIDMapper detects symbol IDs")
        check(mapper.detect_id_type(['ENSG00000141510']) == 'ensembl',
              "GeneIDMapper detects Ensembl IDs")
        check(mapper.detect_id_type(['7157', '672']) == 'entrez',
              "GeneIDMapper detects Entrez IDs")

        result = mapper.convert(['TP53', 'BRCA1'], 'symbol', 'symbol')
        check(result == {'TP53': 'TP53', 'BRCA1': 'BRCA1'},
              "GeneIDMapper same-type conversion (no API)")

        sra_table = pd.DataFrame({
            'sample_alias': ['GSM100001', 'GSM100002', 'SAMN99999'],
            'run_accession': ['SRR1', 'SRR2', 'SRR3'],
        })
        gsm_ids = SRAConnector._extract_gsm_ids(sra_table)
        check(len(gsm_ids) == 2 and 'GSM100001' in gsm_ids,
              "SRAConnector._extract_gsm_ids()")

        ds2 = AcquiredDataset(
            counts_df=acq_counts.rename(columns={c: f'X{c}' for c in acq_counts.columns}),
            metadata=acq_meta.rename(index={c: f'X{c}' for c in acq_meta.index}),
            source_info={'repository': 'test', 'accession': 'TEST002'},
            gene_id_type='symbol',
        )
        engine = PoolingEngine(target_gene_id='symbol', species='Homo sapiens')
        pool = engine.merge([ds, ds2], method='inner')
        check(pool.n_studies == 2 and pool.n_genes == 200,
              f"PoolingEngine.merge() — {pool.n_genes} genes, {pool.n_samples} samples")

        subset = pool.leave_one_study_out(pool.studies[0])
        if isinstance(subset, tuple):
            subset = subset[0]
        check(subset is not None, "PooledDataset.leave_one_study_out()")

        log(f"  [{INFO}] All smoke tests passed")

    except Exception as e:
        check(False, f"Smoke test failed: {e}")
        if verbose:
            log(traceback.format_exc())


# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────

def main():
    start_time = time.time()

    log("\n" + "=" * 65)
    log("  RAPTOR Diagnostic Suite v2.2.2")
    log(f"  Python {sys.version.split()[0]} | {sys.platform}")
    log(f"  Working directory: {Path('.').resolve()}")
    log("=" * 65)

    if not Path('setup.py').exists():
        log(f"\n  [{FAIL}] setup.py not found. Run from the RAPTOR project root.")
        sys.exit(1)

    sys.path.insert(0, str(Path('.').resolve()))

    all_checks = [
        (1,  check_project_structure),
        (2,  check_setup_py),
        (3,  check_imports),
        (4,  check_acquisition_imports),
        (5,  check_class_methods),
        (6,  check_dependencies),
        (7,  check_versions),
        (8,  check_dashboard),
        (9,  check_tests),
        (10, check_cli),
        (11, check_requirements),
        (12, check_code_quality),
        (13, check_git),
        (14, check_cache),
        (15, check_gitignore),
        (16, check_smoke_test),
    ]

    for num, func in all_checks:
        if should_run(num):
            try:
                func()
            except Exception as e:
                log(f"\n  [{FAIL}] Section {num} crashed: {e}")
                if verbose:
                    log(traceback.format_exc())
                results['fail'] += 1

    elapsed = time.time() - start_time
    total = results['pass'] + results['fail'] + results['warn'] + results['skip']

    log(f"\n{'='*65}")
    log(f"  SUMMARY ({elapsed:.1f}s)")
    log(f"{'='*65}")
    log(f"  {results['pass']}/{total} passed  |  "
        f"{results['fail']} failed  |  "
        f"{results['warn']} warnings  |  "
        f"{results['skip']} skipped")

    if results['fail'] == 0:
        log(f"\n  RAPTOR is healthy! Ready to release.\n")
    else:
        log(f"\n  {results['fail']} issue(s) need fixing before release.")
        if not show_fix:
            log(f"     Run with --fix for suggestions.\n")
        else:
            log("")

    if '--save' in sys.argv:
        output_path = Path('check_results.txt')
        clean = re.sub(r'\033\[\d+m', '', all_output.getvalue())
        output_path.write_text(clean, encoding='utf-8')
        print(f"  Results saved to {output_path}")

    return 1 if results['fail'] > 0 else 0


if __name__ == '__main__':
    sys.exit(main())