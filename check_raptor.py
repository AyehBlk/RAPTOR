"""
RAPTOR Diagnostic Suite — Comprehensive health check and debugging tool.

Run routinely before releases, after adding new features/modules, or when
something breaks. Covers structure, imports, dependencies, tests, CLI,
dashboard, code quality, git status, and more.

Usage:
    python check_raptor.py                  # Full diagnostic (all 15 checks)
    python check_raptor.py --quick          # Skip slow checks (tests, CLI, dashboard)
    python check_raptor.py --fix            # Show fix suggestions for every failure
    python check_raptor.py --section 3      # Run only section 3 (imports)
    python check_raptor.py --section 3,5,7  # Run sections 3, 5, and 7
    python check_raptor.py --verbose        # Show extra debug info
    python check_raptor.py --save           # Save results to check_results.txt

How to add a new module:
    1. Add the package path to EXPECTED_PACKAGES
    2. Add module paths to the relevant EXPECTED_*_MODULES list
    3. Add expected exports to the relevant EXPECTED_*_EXPORTS list
    4. Add any new required/optional deps to REQUIRED_DEPS or OPTIONAL_DEPS
    5. Update EXPECTED_TEST_COUNT
    6. Run: python check_raptor.py

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

# -- Packages that must be discoverable by find_packages() or listed in setup.py --
EXPECTED_PACKAGES = [
    'raptor',
    'raptor.utils',
    'raptor.external_modules',
    'raptor.external_modules.acquisition',
    'raptor.dashboard',
    'raptor.dashboard.components',
    'raptor.dashboard.pages',
]

# -- Individual module files that must exist and import --
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

# -- Classes/functions that must be importable from the subpackage __init__ --
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
    'GeneIDMapper': [
        'convert', 'convert_index', 'detect_id_type', 'harmonize_to_common',
    ],
    'PoolingEngine': [
        'merge',
    ],
}

# -- Dependencies --
REQUIRED_DEPS = ['numpy', 'pandas', 'click', 'scipy']

OPTIONAL_DEPS = {
    'pyarrow':    ('Parquet caching (acquisition)', 'pyarrow'),
    'requests':   ('GEO/SRA/TCGA/ArrayExpress API access', 'requests'),
    'GEOparse':   ('GEO dataset parsing', 'GEOparse'),
    'Bio':        ('SRA/GEO Entrez search', 'biopython'),
    'mygene':     ('Gene ID mapping (MyGene.info)', 'mygene'),
    'streamlit':  ('Dashboard', 'streamlit'),
    'plotly':     ('Dashboard visualizations', 'plotly'),
    'sklearn':    ('PCA and normalization', 'scikit-learn'),
}

# -- Dashboard --
EXPECTED_DASHBOARD_PAGES = [
    '01_', '02_', '03_', '04_', '05_',
    '06_', '07_', '08_', '09_', '10_',
]

# -- Tests --
EXPECTED_TEST_COUNT = 105

# -- Version files --
VERSION_LOCATIONS = [
    ('raptor/__init__.py', r'__version__\s*=\s*["\']([^"\']+)["\']'),
    ('setup.py', r'version\s*=\s*["\']([^"\']+)["\']'),
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

    check((root / 'setup.py').exists(), "setup.py exists",
          "You must be in the RAPTOR project root")
    check((root / 'raptor' / '__init__.py').exists(), "raptor/__init__.py exists")
    check((root / 'raptor' / 'external_modules' / 'acquisition' / '__init__.py').exists(),
          "acquisition subpackage __init__.py exists",
          "Create raptor/external_modules/acquisition/__init__.py")
    check((root / 'tests' / 'test_acquisition.py').exists(), "test_acquisition.py exists")

    acq_dir = root / 'raptor' / 'external_modules' / 'acquisition'
    expected_files = [
        'base.py', 'datasets.py', 'cache.py', 'catalog.py',
        'geo.py', 'tcga.py', 'arrayexpress.py', 'sra.py',
        'gene_mapping.py', 'pooling.py',
    ]
    for f in expected_files:
        check((acq_dir / f).exists(), f"acquisition/{f} exists",
              f"Missing file: raptor/external_modules/acquisition/{f}")

    check(
        (root / 'raptor' / 'dashboard' / 'app.py').exists()
        or (root / 'raptor' / 'launch_dashboard.py').exists()
        or (root / 'launch_dashboard.py').exists(),
        "Dashboard entry point exists"
    )


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
                  f"Add '{pkg}' to the packages list in setup.py")

    if 'entry_points' in setup_text or 'console_scripts' in setup_text:
        check(True, "CLI entry point configured")
    else:
        check(False, "CLI entry point (console_scripts) configured",
              "Add entry_points to setup.py", warn_only=True)


# ─────────────────────────────────────────────────────────────────────
# 3. Module Imports
# ─────────────────────────────────────────────────────────────────────

def check_imports():
    section(3, "Module Imports")

    try:
        import raptor
        importlib.reload(raptor)
        check(True, f"import raptor (v{raptor.__version__})")
    except Exception as e:
        check(False, f"import raptor — {e}")
        return

    try:
        from raptor.external_modules import acquisition
        check(True, "from raptor.external_modules import acquisition")
    except Exception as e:
        check(False, f"acquisition import — {e}",
              "Check raptor/external_modules/__init__.py")
        return

    for mod_path in EXPECTED_ACQUISITION_MODULES:
        try:
            importlib.import_module(mod_path)
            check(True, f"import {mod_path.split('.')[-1]}")
        except Exception as e:
            check(False, f"import {mod_path.split('.')[-1]} — {e}")

    missing_exports = []
    for name in EXPECTED_ACQUISITION_EXPORTS:
        try:
            obj = getattr(acquisition, name)
            if obj is None:
                missing_exports.append(f"{name} (is None)")
        except AttributeError:
            missing_exports.append(name)

    check(len(missing_exports) == 0,
          f"All {len(EXPECTED_ACQUISITION_EXPORTS)} expected classes/functions exported",
          f"Missing exports: {missing_exports}" if missing_exports else "")

    # Circular import check
    try:
        for mod_path in EXPECTED_ACQUISITION_MODULES:
            importlib.reload(importlib.import_module(mod_path))
        check(True, "No circular import issues (all modules reloadable)")
    except Exception as e:
        check(False, f"Circular import detected: {e}", warn_only=True)


# ─────────────────────────────────────────────────────────────────────
# 4. Class Methods Verification
# ─────────────────────────────────────────────────────────────────────

def check_class_methods():
    section(4, "Class Methods & Signatures")

    try:
        from raptor.external_modules.acquisition import (
            SRAConnector, GEOConnector, GeneIDMapper, PoolingEngine,
        )
    except ImportError as e:
        check(False, f"Cannot import classes for method check: {e}")
        return

    class_map = {
        'SRAConnector': SRAConnector,
        'GEOConnector': GEOConnector,
        'GeneIDMapper': GeneIDMapper,
        'PoolingEngine': PoolingEngine,
    }

    for cls_name, expected_methods in EXPECTED_CLASS_METHODS.items():
        cls = class_map.get(cls_name)
        if cls is None:
            check(False, f"{cls_name} class not available")
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
# 5. Dependencies
# ─────────────────────────────────────────────────────────────────────

def check_dependencies():
    section(5, "Dependencies")

    for dep in REQUIRED_DEPS:
        try:
            mod = importlib.import_module(dep)
            version = getattr(mod, '__version__', '?')
            check(True, f"{dep} {version} (required)")
        except ImportError:
            check(False, f"{dep} (required) — NOT INSTALLED", f"pip install {dep}")

    for dep, (purpose, pip_name) in OPTIONAL_DEPS.items():
        try:
            mod = importlib.import_module(dep)
            version = getattr(mod, '__version__', '?')
            check(True, f"{dep} {version} — {purpose}")
        except ImportError:
            check(False, f"{dep} — {purpose} — NOT INSTALLED",
                  f"pip install {pip_name}", warn_only=True)


# ─────────────────────────────────────────────────────────────────────
# 6. Version Consistency
# ─────────────────────────────────────────────────────────────────────

def check_versions():
    section(6, "Version Consistency")
    versions = {}

    for filepath, pattern in VERSION_LOCATIONS:
        try:
            text = Path(filepath).read_text(encoding='utf-8')
            match = re.search(pattern, text)
            if match:
                versions[filepath] = match.group(1)
                log(f"  [{INFO}] {filepath}: {match.group(1)}")
            elif filepath == 'setup.py' and '__version__' in text:
                log(f"  [{INFO}] {filepath}: reads version dynamically from __version__ (good practice)")
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

    unique_versions = set(versions.values())
    if len(unique_versions) == 1:
        check(True, f"All versions match: {unique_versions.pop()}")
    elif len(unique_versions) > 1:
        check(False, f"Version MISMATCH: {versions}",
              "Update all files to the same version")
    else:
        check(False, "Could not determine any version")

    # Count hardcoded version references
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
        log(f"  [{INFO}] {hardcoded_count} hardcoded version references found in codebase")


# ─────────────────────────────────────────────────────────────────────
# 7. Dashboard Pages
# ─────────────────────────────────────────────────────────────────────

def check_dashboard():
    section(7, "Dashboard Pages")
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

    # Syntax check all pages
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


# ─────────────────────────────────────────────────────────────────────
# 8. Test Suite
# ─────────────────────────────────────────────────────────────────────

def check_tests():
    section(8, "Test Suite")
    if is_quick():
        log(f"  [{SKIP}] Skipped (--quick mode)")
        results['skip'] += 1
        return

    # Collect
    try:
        result = subprocess.run(
            [sys.executable, '-m', 'pytest', 'tests/test_acquisition.py',
             '--collect-only', '-q'],
            capture_output=True, text=True, timeout=30,
        )
        match = re.search(r'(\d+) test', result.stdout)
        if match:
            n = int(match.group(1))
            check(n >= EXPECTED_TEST_COUNT,
                  f"Collected {n} tests (expected >= {EXPECTED_TEST_COUNT})")
        else:
            check(False, f"Could not parse test count")
    except Exception as e:
        check(False, f"Test collection failed: {e}")

    # Run
    log(f"  [{INFO}] Running full test suite...")
    try:
        result = subprocess.run(
            [sys.executable, '-m', 'pytest', 'tests/test_acquisition.py',
             '-v', '--tb=short', '-q'],
            capture_output=True, text=True, timeout=120,
        )
        passed = int(m.group(1)) if (m := re.search(r'(\d+) passed', result.stdout)) else 0
        failed = int(m.group(1)) if (m := re.search(r'(\d+) failed', result.stdout)) else 0

        check(failed == 0, f"All tests pass ({passed} passed, {failed} failed)",
              "Run: pytest tests/test_acquisition.py -v --tb=long")

        if failed > 0 and verbose:
            for line in result.stdout.split('\n'):
                if 'FAILED' in line:
                    log(f"          {line.strip()}")

        time_match = re.search(r'in ([\d.]+)s', result.stdout)
        if time_match:
            log(f"  [{INFO}] Completed in {time_match.group(1)}s")

    except subprocess.TimeoutExpired:
        check(False, "Test run timed out (>120s)")
    except Exception as e:
        check(False, f"Test run failed: {e}")


# ─────────────────────────────────────────────────────────────────────
# 9. CLI Entry Point
# ─────────────────────────────────────────────────────────────────────

def check_cli():
    section(9, "CLI Entry Point")
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

            # Compare CLI version to __init__ version
            try:
                import raptor
                init_ver = raptor.__version__
                if init_ver not in cli_output:
                    check(False,
                          f"CLI version mismatch: CLI says '{cli_output}' but __init__ is '{init_ver}'",
                          "Update cli.py to use: @click.version_option(version=raptor.__version__)",
                          warn_only=True)
            except Exception:
                pass
        else:
            check(False, "raptor --version failed",
                  "Check CLI entry point in setup.py", warn_only=True)
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
# 10. Requirements & Environment Files
# ─────────────────────────────────────────────────────────────────────

def check_requirements():
    section(10, "Requirements & Environment Files")
    req_file = Path('requirements.txt')
    if not req_file.exists():
        check(False, "requirements.txt exists")
        return

    req_text = req_file.read_text(encoding='utf-8').lower()
    for dep in ['numpy', 'pandas', 'click', 'scipy']:
        check(dep in req_text, f"requirements.txt includes {dep}")

    for dep in ['requests', 'pyarrow']:
        check(dep in req_text, f"requirements.txt includes {dep} (acquisition)",
              f"Add {dep} to requirements.txt", warn_only=True)

    for env_file in ['environment.yml', 'environment-full.yml']:
        if Path(env_file).exists():
            log(f"  [{INFO}] {env_file} exists ({Path(env_file).stat().st_size:,} bytes)")


# ─────────────────────────────────────────────────────────────────────
# 11. Code Quality
# ─────────────────────────────────────────────────────────────────────

def check_code_quality():
    section(11, "Code Quality")

    # UTF-8 encoding check
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

    # print() in library code
    print_locs = []
    for pyfile in Path('raptor/external_modules').rglob('*.py'):
        try:
            for i, line in enumerate(pyfile.read_text(encoding='utf-8').split('\n'), 1):
                if line.strip().startswith('print(') and not line.strip().startswith('#'):
                    print_locs.append(f"{pyfile.name}:{i}")
        except Exception:
            pass

    check(len(print_locs) == 0, "No print() in library code (use logging instead)",
          f"Found in: {print_locs[:5]}", warn_only=True)

    # Syntax check all files
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

    # Missing encoding='utf-8' in file writes
    missing_enc = []
    for pyfile in Path('raptor/external_modules').rglob('*.py'):
        try:
            content = pyfile.read_text(encoding='utf-8')
            for i, line in enumerate(content.split('\n'), 1):
                if ('open(' in line and ("'w'" in line or '"w"' in line)
                        and 'encoding' not in line and not line.strip().startswith('#')):
                    missing_enc.append(f"{pyfile.name}:{i}")
        except Exception:
            pass

    if missing_enc:
        check(False, "All file writes use encoding='utf-8'",
              f"Missing in: {missing_enc[:5]}", warn_only=True)
    else:
        check(True, "No obvious missing encoding='utf-8' in write operations")


# ─────────────────────────────────────────────────────────────────────
# 12. Git Repository
# ─────────────────────────────────────────────────────────────────────

def check_git():
    section(12, "Git Repository")
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
# 13. Cache & Data Directories
# ─────────────────────────────────────────────────────────────────────

def check_cache():
    section(13, "Cache & Data Directories")
    cache_dir = Path.home() / '.raptor_cache'
    if cache_dir.exists():
        n_files = sum(1 for _ in cache_dir.rglob('*') if _.is_file())
        total_size = sum(f.stat().st_size for f in cache_dir.rglob('*') if f.is_file())
        log(f"  [{INFO}] Cache: {cache_dir} ({n_files} files, {total_size / 1024 / 1024:.1f} MB)")
    else:
        log(f"  [{INFO}] No cache directory yet (created on first download)")

    pycache_dirs = list(Path('raptor').rglob('__pycache__'))
    log(f"  [{INFO}] {len(pycache_dirs)} __pycache__ directories")
    check(True, "Cache check complete")


# ─────────────────────────────────────────────────────────────────────
# 14. .gitignore
# ─────────────────────────────────────────────────────────────────────

def check_gitignore():
    section(14, ".gitignore")
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
# 15. Smoke Test (offline functionality)
# ─────────────────────────────────────────────────────────────────────

def check_smoke_test():
    section(15, "Smoke Test (offline functionality)")
    if is_quick():
        log(f"  [{SKIP}] Skipped (--quick mode)")
        results['skip'] += 1
        return

    try:
        import numpy as np
        import pandas as pd
        from raptor.external_modules.acquisition import (
            AcquiredDataset, PooledDataset, CacheManager, GeneIDMapper,
            PoolingEngine, SRAConnector,
        )

        # Create a mini dataset (200 genes — above PoolingEngine min threshold)
        counts = pd.DataFrame(
            np.random.randint(0, 100, (200, 4)),
            index=[f'GENE_{i}' for i in range(200)],
            columns=['S1', 'S2', 'S3', 'S4'],
        )
        counts.index.name = 'gene_id'
        meta = pd.DataFrame(
            {'condition': ['ctrl', 'ctrl', 'treat', 'treat']},
            index=['S1', 'S2', 'S3', 'S4'],
        )
        meta.index.name = 'sample_id'

        ds = AcquiredDataset(
            counts_df=counts, metadata=meta,
            source_info={'repository': 'test', 'accession': 'TEST001'},
            gene_id_type='symbol',
        )
        check(ds.n_genes == 200 and ds.n_samples == 4, "AcquiredDataset creation")

        # Integrity validation
        integrity = ds.validate_integrity()
        check(integrity['valid'], "AcquiredDataset.validate_integrity()")

        # GeneIDMapper detection
        mapper = GeneIDMapper('Homo sapiens')
        check(mapper.detect_id_type(['TP53', 'BRCA1', 'EGFR']) == 'symbol',
              "GeneIDMapper.detect_id_type()")
        check(mapper.detect_id_type(['ENSG00000141510']) == 'ensembl',
              "GeneIDMapper detects Ensembl IDs")
        check(mapper.detect_id_type(['7157', '672']) == 'entrez',
              "GeneIDMapper detects Entrez IDs")

        # Same-type conversion (no API needed)
        result = mapper.convert(['TP53', 'BRCA1'], 'symbol', 'symbol')
        check(result == {'TP53': 'TP53', 'BRCA1': 'BRCA1'},
              "GeneIDMapper same-type conversion (no API)")

        # SRA GSM extraction
        sra_table = pd.DataFrame({
            'sample_alias': ['GSM100001', 'GSM100002', 'SAMN99999'],
            'run_accession': ['SRR1', 'SRR2', 'SRR3'],
        })
        gsm_ids = SRAConnector._extract_gsm_ids(sra_table)
        check(len(gsm_ids) == 2 and 'GSM100001' in gsm_ids,
              "SRAConnector._extract_gsm_ids()")

        # Pooling two datasets
        ds2 = AcquiredDataset(
            counts_df=counts.rename(columns={c: f'X{c}' for c in counts.columns}),
            metadata=meta.rename(index={c: f'X{c}' for c in meta.index}),
            source_info={'repository': 'test', 'accession': 'TEST002'},
            gene_id_type='symbol',
        )
        engine = PoolingEngine(target_gene_id='symbol', species='Homo sapiens')
        pool = engine.merge([ds, ds2], method='inner')
        check(pool.n_studies == 2 and pool.n_genes == 200,
              f"PoolingEngine.merge() — {pool.n_genes} genes, {pool.n_samples} samples")

        # Leave-one-study-out
        subset = pool.leave_one_study_out(pool.studies[0])
        # May return PooledDataset or tuple — handle both
        if isinstance(subset, tuple):
            subset = subset[0]
        check(subset is not None,
              "PooledDataset.leave_one_study_out()")

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
    log("  RAPTOR Diagnostic Suite")
    log(f"  Python {sys.version.split()[0]} | {sys.platform}")
    log(f"  Working directory: {Path('.').resolve()}")
    log("=" * 65)

    if not Path('setup.py').exists():
        log(f"\n  [{FAIL}] setup.py not found. Run from the RAPTOR project root.")
        sys.exit(1)

    sys.path.insert(0, str(Path('.').resolve()))

    all_checks = [
        (1, check_project_structure),   (2, check_setup_py),
        (3, check_imports),             (4, check_class_methods),
        (5, check_dependencies),        (6, check_versions),
        (7, check_dashboard),           (8, check_tests),
        (9, check_cli),                 (10, check_requirements),
        (11, check_code_quality),       (12, check_git),
        (13, check_cache),              (14, check_gitignore),
        (15, check_smoke_test),
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
