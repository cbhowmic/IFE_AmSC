import { memo, startTransition, useCallback, useDeferredValue, useEffect, useMemo, useState } from 'react';
import {
  Background,
  BackgroundVariant,
  Controls,
  Handle,
  MarkerType,
  MiniMap,
  Panel,
  Position,
  ReactFlow,
  applyNodeChanges,
} from '@xyflow/react';
import '@xyflow/react/dist/style.css';
import './App.css';

const PLAYBACK_OPTIONS = [
  { label: 'Slow', frameStep: 4, intervalMs: 110 },
  { label: 'Medium', frameStep: 10, intervalMs: 80 },
  { label: 'Fast', frameStep: 24, intervalMs: 45 },
];

const DATA_FILES = {
  payload: 'graph_payload.json',
  nodes: 'nodes_rf.json',
  edges: 'edges_rf.json',
  time: 'time_rf.json',
};

const DATA_TIMEOUT_MS = 5500;
const CANVAS_BASE = '#182533';
const CANVAS_GRID = 'rgba(129, 169, 203, 0.1)';
const SURFACE = 'rgba(13, 26, 38, 0.92)';
const SURFACE_BORDER = 'rgba(120, 171, 218, 0.26)';
const FLOW_COLOR = '#5a9fe4';
const LOSS_COLOR = '#242a31';
const LOSS_EDGE_COLOR = '#f0aa47';
const UNSPECIFIED_COLOR = '#dd8a43';
const LOSS_NODE_COLOR = LOSS_COLOR;
const NODE_FACE = '#e7f1f7';
const NODE_FACE_LOSS = '#cfd5d8';
const NODE_FACE_DIAGNOSTIC = '#f5dfca';
const NODE_FRAME = '#173651';
const PANEL_TEXT = '#dcecff';
const WARNING_AMBER = '#f0aa47';
const WARNING_RED = '#df684b';
const WARNING_PURPLE = '#7b68b6';
const HANDLE_COLOR = '#050607';
const HANDLE_BORDER = '#dde3e8';
const ALERT_WINDOW_FRAMES = 28;

const ALERT_PRESENTATION = {
  none: { color: NODE_FRAME, label: '' },
  half: { color: WARNING_AMBER, label: '0.5 steady state' },
  double: { color: WARNING_RED, label: '2x initial mass' },
  both: { color: WARNING_PURPLE, label: 'dual threshold' },
};

const SURROGATE_DEFAULTS = {
  isotopeSeparationTritium: {
    label: 'Isotope separation T',
    value: null,
    unit: 'g',
  },
  plantDoublingTime: {
    label: 'Plant doubling time',
    value: null,
    unit: 'd',
  },
  minimalStartupInventory: {
    label: 'Minimal startup inventory',
    value: null,
    unit: 'g',
  },
};

function normalizeBaseUrl(value) {
  return typeof value === 'string' ? value.trim().replace(/\/+$/, '') : '';
}

function getRemoteBaseUrl() {
  const queryBase = new URLSearchParams(window.location.search).get('dataBaseUrl');
  return normalizeBaseUrl(queryBase || import.meta.env.VITE_DATA_BASE_URL || '');
}

function buildDataUrl(baseUrl, fileName) {
  return baseUrl ? `${baseUrl}/${fileName}` : `/${fileName}`;
}

async function fetchJsonWithTimeout(url, timeoutMs = DATA_TIMEOUT_MS) {
  const controller = new AbortController();
  const timeout = window.setTimeout(() => controller.abort(), timeoutMs);

  try {
    const response = await fetch(url, { signal: controller.signal });

    if (!response.ok) {
      throw new Error(`${response.status} ${response.statusText}`);
    }

    return response.json();
  } finally {
    window.clearTimeout(timeout);
  }
}

async function loadJsonWithFallback(fileName, remoteBaseUrl) {
  if (remoteBaseUrl) {
    try {
      return {
        data: await fetchJsonWithTimeout(buildDataUrl(remoteBaseUrl, fileName)),
        source: 'remote',
        error: null,
      };
    } catch (error) {
      console.warn(`Remote ${fileName} failed; falling back to bundled JSON.`, error);
    }
  }

  return {
    data: await fetchJsonWithTimeout(buildDataUrl('', fileName)),
    source: remoteBaseUrl ? 'fallback' : 'local',
    error: null,
  };
}

function getDatasetSource(results, remoteBaseUrl, sourceKind = 'JSON') {
  if (!remoteBaseUrl) {
    return `Bundled ${sourceKind}`;
  }

  if (results.every((result) => result.source === 'remote')) {
    return `Remote ${sourceKind}`;
  }

  if (results.every((result) => result.source === 'fallback')) {
    return 'Bundled fallback';
  }

  return 'Mixed remote/local';
}

async function loadVisualizationSources(remoteBaseUrl) {
  try {
    const payloadResult = await loadJsonWithFallback(DATA_FILES.payload, remoteBaseUrl);
    const payload = payloadResult.data;

    if (
      payload &&
      typeof payload === 'object' &&
      Array.isArray(payload.nodes) &&
      Array.isArray(payload.edges)
    ) {
      return {
        nodesResult: { data: payload.nodes, source: payloadResult.source },
        edgesResult: { data: payload.edges, source: payloadResult.source },
        timeResult: { data: payload, source: payloadResult.source },
        results: [payloadResult],
        sourceKind: 'payload',
      };
    }

    throw new Error(`${DATA_FILES.payload} is missing nodes or edges`);
  } catch (error) {
    console.warn(`${DATA_FILES.payload} unavailable; falling back to split JSON files.`, error);
  }

  const [nodesResult, edgesResult, timeResult] = await Promise.all([
    loadJsonWithFallback(DATA_FILES.nodes, remoteBaseUrl),
    loadJsonWithFallback(DATA_FILES.edges, remoteBaseUrl),
    loadJsonWithFallback(DATA_FILES.time, remoteBaseUrl),
  ]);

  return {
    nodesResult,
    edgesResult,
    timeResult,
    results: [nodesResult, edgesResult, timeResult],
    sourceKind: 'JSON',
  };
}

function findFirstThresholdCrossing(series, threshold) {
  if (!Array.isArray(series) || !Number.isFinite(threshold)) {
    return -1;
  }

  for (let index = 1; index < series.length; index += 1) {
    const previous = series[index - 1];
    const current = series[index];

    if (!Number.isFinite(previous) || !Number.isFinite(current)) {
      continue;
    }

    if (previous < threshold && current >= threshold) {
      return index;
    }
  }

  return -1;
}

function getAlertBeepLevel(currentIndex, crossingIndex) {
  if (!Number.isFinite(currentIndex) || crossingIndex < 0) {
    return 0;
  }

  const elapsed = currentIndex - crossingIndex;
  if (elapsed < 0 || elapsed > ALERT_WINDOW_FRAMES) {
    return 0;
  }

  if ((elapsed >= 0 && elapsed <= 4) || (elapsed >= 10 && elapsed <= 14)) {
    return 2;
  }

  if ((elapsed >= 5 && elapsed <= 7) || (elapsed >= 15 && elapsed <= 18)) {
    return 1;
  }

  return 0;
}

function formatMass(value) {
  if (!Number.isFinite(value)) {
    return 'n/a';
  }

  if (Math.abs(value) >= 1000) {
    return value.toLocaleString(undefined, { maximumFractionDigits: 0 });
  }

  if (Math.abs(value) >= 10) {
    return value.toLocaleString(undefined, { maximumFractionDigits: 1 });
  }

  return value.toLocaleString(undefined, { maximumFractionDigits: 3 });
}

function formatPercent(value) {
  if (!Number.isFinite(value)) {
    return 'n/a';
  }

  return `${Math.round(value * 100)}%`;
}

function formatPredictionValue(prediction) {
  const value = Number(prediction?.value);

  if (!Number.isFinite(value)) {
    return 'pending';
  }

  const formatted = formatMass(value);
  if (!prediction?.unit) {
    return formatted;
  }

  const unit = String(prediction.unit).trim();
  return unit.startsWith('[') ? `${formatted} ${unit}` : `${formatted} [${unit}]`;
}

function formatMassWithUnit(value, unit = 'g') {
  const formatted = formatMass(value);
  return formatted === 'n/a' ? formatted : `${formatted} [${unit}]`;
}

function formatTimeValue(value, unit) {
  if (!Number.isFinite(value)) {
    return 'Time n/a';
  }

  const rounded =
    Math.abs(value) >= 100 ? value.toLocaleString(undefined, { maximumFractionDigits: 0 }) : value.toLocaleString(undefined, { maximumFractionDigits: 2 });

  return `Time ${rounded} ${unit}`;
}

function findClosestTimeIndex(series, targetValue) {
  if (!Array.isArray(series) || !series.length || !Number.isFinite(targetValue)) {
    return -1;
  }

  let bestIndex = 0;
  let bestDistance = Math.abs(Number(series[0]) - targetValue);

  for (let index = 1; index < series.length; index += 1) {
    const value = Number(series[index]);

    if (!Number.isFinite(value)) {
      continue;
    }

    const distance = Math.abs(value - targetValue);
    if (distance < bestDistance) {
      bestDistance = distance;
      bestIndex = index;
    }
  }

  return bestIndex;
}

function getMassColor(fillFrac, isLossNode) {
  if (isLossNode) {
    return '#5f5b66';
  }

  if (fillFrac < 0.25) {
    return '#e59b48';
  }

  if (fillFrac < 0.5) {
    return '#8fb8dd';
  }

  if (fillFrac < 0.85) {
    return '#4f91d2';
  }

  return '#2364af';
}

function isDiagnosticLabel(label) {
  return typeof label === 'string' && label.toLowerCase().includes('box');
}

function formatSubsystemLabel(label) {
  return String(label || '').replaceAll('_', ' ');
}

function isIsotopeSeparationLabel(label) {
  const normalized = String(label || '').toLowerCase().replace(/[^a-z]/g, '');
  return normalized.includes('isotopeseparation') || normalized.includes('isotopeseperation');
}

function normalizeSurrogatePredictions(value) {
  const source = value && typeof value === 'object' ? value : {};

  return Object.fromEntries(
    Object.entries(SURROGATE_DEFAULTS).map(([key, fallback]) => {
      const raw = source[key];

      if (typeof raw === 'number') {
        return [key, { ...fallback, value: raw }];
      }

      if (raw && typeof raw === 'object') {
        return [key, { ...fallback, ...raw }];
      }

      return [key, fallback];
    }),
  );
}

function normalizeTimePayload(rawTimeData) {
  if (!rawTimeData || typeof rawTimeData !== 'object') {
    return {
      timeSeries: [],
      timeUnit: 'days',
      surrogatePredictions: SURROGATE_DEFAULTS,
    };
  }

  const timeData = rawTimeData.time && typeof rawTimeData.time === 'object'
    ? rawTimeData.time
    : rawTimeData;

  return {
    timeSeries: Array.isArray(timeData.timeSeries) ? timeData.timeSeries : [],
    timeUnit: typeof timeData.timeUnit === 'string' ? timeData.timeUnit : 'days',
    surrogatePredictions: normalizeSurrogatePredictions(
      rawTimeData.surrogatePredictions ?? timeData.surrogatePredictions,
    ),
  };
}

function getAlertKind(hasDoubleInitial, hasHalfSteady) {
  if (hasDoubleInitial && hasHalfSteady) {
    return 'both';
  }

  if (hasDoubleInitial) {
    return 'double';
  }

  if (hasHalfSteady) {
    return 'half';
  }

  return 'none';
}

function getEdgePresentation(edge) {
  if (edge.data?.isLoss) {
    return {
      stroke: LOSS_EDGE_COLOR,
      strokeWidth: 4.2,
      opacity: 0.95,
      animated: false,
    };
  }

  if (edge.data?.isUnknownFlow) {
    return {
      stroke: UNSPECIFIED_COLOR,
      strokeWidth: 2.2,
      opacity: 0.72,
      animated: false,
    };
  }

  const baseFlow = Number(edge.data?.baseFlow);
  const normalizedFlow = Number.isFinite(baseFlow)
    ? Math.max(0, Math.min(1, baseFlow))
    : 0.5;

  return {
    stroke: FLOW_COLOR,
    strokeWidth: 1.8 + normalizedFlow * 3.8,
    opacity: 0.45 + normalizedFlow * 0.55,
    animated: true,
  };
}

function getBestHandlePair(sourcePosition, targetPosition, isLossEdge = false) {
  if (!sourcePosition || !targetPosition) {
    return { sourceHandle: 'right', targetHandle: 'left' };
  }

  if (isLossEdge) {
    return targetPosition.y >= sourcePosition.y
      ? { sourceHandle: 'bottom', targetHandle: 'top' }
      : { sourceHandle: 'top', targetHandle: 'bottom' };
  }

  const dx = targetPosition.x - sourcePosition.x;
  const dy = targetPosition.y - sourcePosition.y;

  if (Math.abs(dx) >= Math.abs(dy)) {
    return dx >= 0
      ? { sourceHandle: 'right', targetHandle: 'left' }
      : { sourceHandle: 'left', targetHandle: 'right' };
  }

  return dy >= 0
    ? { sourceHandle: 'bottom', targetHandle: 'top' }
    : { sourceHandle: 'top', targetHandle: 'bottom' };
}

const TARGET_HANDLES = [
  { id: 'left', position: Position.Left },
  { id: 'top', position: Position.Top },
  { id: 'bottom', position: Position.Bottom },
];

const SOURCE_HANDLES = [
  { id: 'right', position: Position.Right },
  { id: 'left', position: Position.Left },
  { id: 'top', position: Position.Top },
  { id: 'bottom', position: Position.Bottom },
];

function handleStyle(background) {
  return {
    background,
    width: 10,
    height: 10,
    border: `2px solid ${HANDLE_BORDER}`,
    borderRadius: 0,
  };
}

const ComponentNode = memo(function ComponentNode({ data }) {
  const {
    alertKind,
    alertPulseLevel,
    alertText,
    currentMass,
    currentToSteadyRatio,
    fillFrac,
    hasMeaningfulSteadyState,
    initialMass,
    isDiagnostic,
    isLossNode,
    isSurrogateTarget,
    label,
    maxMass,
    minMass,
    steadyMass,
  } = data;

  const presentation = ALERT_PRESENTATION[alertKind] ?? ALERT_PRESENTATION.none;
  const hasAlert = !isLossNode && !isDiagnostic && alertKind !== 'none';
  const alertClass = hasAlert ? `is-${alertKind}` : 'is-calm';
  const pulseClass = alertPulseLevel > 1 ? 'is-pulsing-strong' : alertPulseLevel > 0 ? 'is-pulsing' : '';
  const nodeFaceColor = isLossNode ? NODE_FACE_LOSS : isDiagnostic ? NODE_FACE_DIAGNOSTIC : NODE_FACE;
  const frameColor = isLossNode ? LOSS_COLOR : isDiagnostic ? UNSPECIFIED_COLOR : presentation.color;
  const secondStatLabel = isDiagnostic ? 'Readout span' : 'Fullness';

  return (
    <div
      className={`plant-node-shell ${alertClass} ${pulseClass} ${isDiagnostic ? 'is-diagnostic' : ''} ${isSurrogateTarget ? 'is-surrogate-target' : ''}`}
      style={{
        '--node-fill-height': `${Math.max(4, fillFrac * 100)}%`,
        '--node-fill-color': getMassColor(fillFrac, isLossNode),
        '--node-face-color': nodeFaceColor,
        '--node-fullness-color': getMassColor(fillFrac, isLossNode),
        '--node-frame-color': frameColor,
        '--node-alert-color': presentation.color,
      }}
    >
      <div className={`plant-node-card ${isLossNode ? 'is-loss-node' : ''} ${isDiagnostic ? 'is-diagnostic-node' : ''}`}>
        {hasAlert && <div className="plant-node-alert-chip" />}

        {TARGET_HANDLES.map((handle) => (
          <Handle
            key={`target-${handle.id}`}
            id={handle.id}
            type="target"
            position={handle.position}
            style={handleStyle(HANDLE_COLOR)}
          />
        ))}
        {!isLossNode && (
          <>
            {SOURCE_HANDLES.map((handle) => (
              <Handle
                key={`source-${handle.id}`}
                id={handle.id}
                type="source"
                position={handle.position}
                style={handleStyle(FLOW_COLOR)}
              />
            ))}
          </>
        )}

        <div className="plant-node-fill" />

        <div className="plant-node-content">
          <div className="plant-node-title">
            {label}
          </div>
          {isSurrogateTarget && (
            <div className="surrogate-target-tag">Surrogate target</div>
          )}

          <div className="plant-node-stats">
            <div className="plant-node-stat">
              <div className="plant-node-stat-label">
              Current mass
              </div>
              <div className="plant-node-stat-value">
                {formatMass(currentMass)}
              </div>
            </div>

            <div className="plant-node-stat">
              <div className="plant-node-stat-label">
              {secondStatLabel}
              </div>
              <div className="plant-node-stat-value">
                {formatPercent(fillFrac)}
              </div>
            </div>
          </div>

          <div className="plant-node-foot">
            <span>init {formatMass(initialMass)}</span>
            {hasMeaningfulSteadyState ? (
              <>
                <span>ss {formatMass(steadyMass)}</span>
                <span>{formatPercent(currentToSteadyRatio)}</span>
              </>
            ) : (
              <span>diagnostic</span>
            )}
          </div>

          <div className="plant-node-range">
            {isDiagnostic ? 'signal' : 'range'} {formatMass(minMass)} to {formatMass(maxMass)}
          </div>

          <div className="plant-node-badges">
            {alertText && (
              <span className="plant-node-badge">
                {alertText}
              </span>
              )}
          </div>
        </div>
      </div>
    </div>
  );
});

const nodeTypes = {
  componentNode: ComponentNode,
};

export default function App() {
  const [baseNodes, setBaseNodes] = useState([]);
  const [allEdges, setAllEdges] = useState([]);
  const [massLookup, setMassLookup] = useState({});
  const [timeSeries, setTimeSeries] = useState([]);
  const [timeUnit, setTimeUnit] = useState('days');
  const [surrogatePredictions, setSurrogatePredictions] = useState(SURROGATE_DEFAULTS);
  const [timeIndex, setTimeIndex] = useState(0);
  const [isPlaying, setIsPlaying] = useState(false);
  const [playbackMode, setPlaybackMode] = useState(1);
  const [layoutStatus, setLayoutStatus] = useState('');
  const [showLossEdges, setShowLossEdges] = useState(true);
  const [showUnspecifiedEdges, setShowUnspecifiedEdges] = useState(false);
  const [dataStatus, setDataStatus] = useState(() => {
    const remoteBaseUrl = getRemoteBaseUrl();

    return {
      state: 'loading',
      source: remoteBaseUrl ? 'Remote JSON' : 'Bundled JSON',
      remoteBaseUrl,
      message: remoteBaseUrl ? `Trying ${remoteBaseUrl}` : 'Loading bundled JSON',
    };
  });

  const loadVisualizationData = useCallback(async (forceBundled = false, shouldApply = () => true) => {
    const configuredRemoteBaseUrl = getRemoteBaseUrl();
    const remoteBaseUrl = forceBundled ? '' : configuredRemoteBaseUrl;

    setDataStatus({
      state: 'loading',
      source: forceBundled ? 'Bundled override' : remoteBaseUrl ? 'Remote JSON' : 'Bundled JSON',
      remoteBaseUrl: configuredRemoteBaseUrl,
      message: forceBundled
        ? 'Loading bundled fallback files'
        : remoteBaseUrl
          ? `Trying ${remoteBaseUrl}`
          : 'Loading bundled JSON',
    });

    try {
      const {
        nodesResult,
        edgesResult,
        timeResult,
        results,
        sourceKind,
      } = await loadVisualizationSources(remoteBaseUrl);

      if (!shouldApply()) {
        return;
      }

      const nodesData = Array.isArray(nodesResult.data) ? nodesResult.data : [];
      const edgesData = Array.isArray(edgesResult.data) ? edgesResult.data : [];
      const timeData = normalizeTimePayload(timeResult.data);
      const lookup = {};
      const positionLookup = {};
      const typedNodes = nodesData.map((node) => {
        const nodeData = node.data || {};
        const rawLabel = nodeData.label || `Subsystem ${node.id}`;
        const label = formatSubsystemLabel(rawLabel);
        const isDiagnostic = Boolean(nodeData.isDiagnostic) || isDiagnosticLabel(label);
        const isSurrogateTarget = isIsotopeSeparationLabel(rawLabel) || isIsotopeSeparationLabel(label);
        const hasMeaningfulSteadyState = Boolean(
          nodeData.hasMeaningfulSteadyState ?? !isDiagnostic,
        );
        const series = Array.isArray(nodeData.massSeries) ? nodeData.massSeries : [];
        const initialMass = Number(nodeData.initialMass) || 0;
        const finiteSeries = series.filter(Number.isFinite);
        const steadyMassFromSeries = finiteSeries.at(-1);
        const steadyMass = Number.isFinite(nodeData.steadyMass)
          ? nodeData.steadyMass
          : steadyMassFromSeries ?? initialMass;
        const halfSteadyThreshold =
          hasMeaningfulSteadyState && steadyMass > 0 ? 0.5 * steadyMass : Number.NaN;
        const doubleInitialThreshold =
          !isDiagnostic && initialMass > 0 ? 2 * initialMass : Number.NaN;

        lookup[node.id] = series;
        positionLookup[node.id] = node.position;

        return {
          ...node,
          type: 'componentNode',
          data: {
            label,
            initialMass,
            minMass: finiteSeries.length ? Math.min(...finiteSeries) : initialMass,
            maxMass: finiteSeries.length ? Math.max(...finiteSeries) : initialMass,
            steadyMass,
            halfSteadyCrossingIndex: findFirstThresholdCrossing(series, halfSteadyThreshold),
            doubleInitialCrossingIndex: findFirstThresholdCrossing(series, doubleInitialThreshold),
            currentMass: initialMass,
            currentToSteadyRatio: steadyMass > 0 ? initialMass / steadyMass : 0,
            fillFrac: 0,
            alertKind: 'none',
            alertPulseLevel: 0,
            alertText: '',
            isLossNode: Boolean(nodeData.isLossNode),
            isDiagnostic,
            isSurrogateTarget,
            hasMeaningfulSteadyState,
          },
        };
      });

      const styledEdges = edgesData.map((edge) => {
        const presentation = getEdgePresentation(edge);
        const handles = getBestHandlePair(
          positionLookup[edge.source],
          positionLookup[edge.target],
          edge.data?.isLoss,
        );

        return {
          ...edge,
          type: 'smoothstep',
          sourceHandle: handles.sourceHandle,
          targetHandle: handles.targetHandle,
          animated: presentation.animated,
          style: {
            stroke: presentation.stroke,
            strokeWidth: presentation.strokeWidth,
            opacity: presentation.opacity,
          },
          markerEnd: {
            type: MarkerType.ArrowClosed,
            width: 18,
            height: 18,
            color: presentation.stroke,
          },
        };
      });

      setMassLookup(lookup);
      setBaseNodes(typedNodes);
      setAllEdges(styledEdges);
      setTimeSeries(timeData.timeSeries);
      setTimeUnit(timeData.timeUnit);
      setSurrogatePredictions(timeData.surrogatePredictions);
      setTimeIndex(0);
      setIsPlaying(false);
      setDataStatus({
        state: 'ready',
        source: forceBundled
          ? 'Bundled override'
          : getDatasetSource(results, remoteBaseUrl, sourceKind),
        remoteBaseUrl: configuredRemoteBaseUrl,
        message: `${typedNodes.length} nodes, ${styledEdges.length} edges, ${timeData.timeSeries.length} time samples`,
      });
    } catch (error) {
      if (!shouldApply()) {
        return;
      }

      console.error('Unable to load visualization data.', error);
      setDataStatus({
        state: 'error',
        source: 'Data unavailable',
        remoteBaseUrl: configuredRemoteBaseUrl,
        message: error instanceof Error ? error.message : 'Unknown data loading error',
      });
    }
  }, []);

  useEffect(() => {
    let isMounted = true;
    window.queueMicrotask(() => {
      loadVisualizationData(false, () => isMounted);
    });

    return () => {
      isMounted = false;
    };
  }, [loadVisualizationData]);

  const deferredTimeIndex = useDeferredValue(timeIndex);
  const playback = PLAYBACK_OPTIONS[playbackMode];

  const maxTime = useMemo(() => {
    const firstNode = baseNodes[0];
    const firstSeries = firstNode ? massLookup[firstNode.id] : null;
    return firstSeries?.length ? firstSeries.length - 1 : 0;
  }, [baseNodes, massLookup]);

  const boundedTimeIndex = Math.min(timeIndex, maxTime);
  const currentTimeValue = timeSeries[boundedTimeIndex] ?? boundedTimeIndex;
  const maxTimeValue = timeSeries[maxTime] ?? maxTime;

  useEffect(() => {
    if (!isPlaying || maxTime <= 0) {
      return undefined;
    }

    const timer = setInterval(() => {
      startTransition(() => {
        setTimeIndex((current) => {
          const next = current + playback.frameStep;
          return next >= maxTime ? 0 : next;
        });
      });
    }, playback.intervalMs);

    return () => clearInterval(timer);
  }, [isPlaying, maxTime, playback]);

  const deferredBoundedTimeIndex = Math.min(deferredTimeIndex, maxTime);

  const nodes = useMemo(() => {
    return baseNodes.map((node) => {
      const series = massLookup[node.id] || [];
      const currentMass =
        series[Math.min(deferredBoundedTimeIndex, Math.max(series.length - 1, 0))] ??
        node.data.initialMass;
      const steadyMass = node.data.steadyMass;
      const minMass = node.data.minMass;
      const maxMass = node.data.maxMass;
      const halfSteadyAlert = getAlertBeepLevel(
        deferredBoundedTimeIndex,
        node.data.halfSteadyCrossingIndex,
      );
      const doubleInitialAlert = getAlertBeepLevel(
        deferredBoundedTimeIndex,
        node.data.doubleInitialCrossingIndex,
      );
      const alertPulseLevel = Math.max(halfSteadyAlert, doubleInitialAlert);
      const hasHalfSteady =
        !node.data.isLossNode &&
        !node.data.isDiagnostic &&
        node.data.hasMeaningfulSteadyState &&
        Number.isFinite(steadyMass) &&
        steadyMass > 0 &&
        currentMass >= 0.5 * steadyMass;
      const hasDoubleInitial =
        !node.data.isLossNode &&
        !node.data.isDiagnostic &&
        Number.isFinite(node.data.initialMass) &&
        node.data.initialMass > 0 &&
        currentMass >= 2 * node.data.initialMass;
      const alertKind = getAlertKind(hasDoubleInitial, hasHalfSteady);
      let alertText = '';

      if (alertKind === 'both') {
        alertText = 'dual threshold';
      } else if (alertKind === 'double') {
        alertText = 'reached 2x init';
      } else if (alertKind === 'half') {
        alertText = 'reached 0.5 ss';
      }

      const currentToSteadyRatio =
        node.data.hasMeaningfulSteadyState && Number.isFinite(steadyMass) && steadyMass > 0
          ? currentMass / steadyMass
          : 0;
      const fillRange = Number.isFinite(maxMass) && Number.isFinite(minMass) ? maxMass - minMass : 0;
      const fillFrac = fillRange > 0
        ? Math.max(0, Math.min(1, (currentMass - minMass) / fillRange))
        : 0;

      return {
        ...node,
        data: {
          ...node.data,
          alertKind,
          alertPulseLevel,
          alertText,
          currentMass,
          currentToSteadyRatio,
          fillFrac,
        },
      };
    });
  }, [baseNodes, deferredBoundedTimeIndex, massLookup]);

  const onNodesChange = useCallback((changes) => {
    setBaseNodes((existingNodes) => applyNodeChanges(changes, existingNodes));
  }, []);

  const massStats = useMemo(() => {
    let totalMass = 0;
    let steadyStateMass = 0;
    let highestRatio = 0;
    let aboveDouble = 0;
    let aboveHalfSteady = 0;

    for (const node of nodes) {
      if (node.data.isLossNode || node.data.isDiagnostic) {
        continue;
      }

      totalMass += Number(node.data.currentMass) || 0;
      if (node.data.hasMeaningfulSteadyState) {
        steadyStateMass += Number(node.data.steadyMass) || 0;
        highestRatio = Math.max(highestRatio, node.data.currentToSteadyRatio || 0);
      }

      if (
        Number.isFinite(node.data.initialMass) &&
        node.data.initialMass > 0 &&
        node.data.currentMass >= 2 * node.data.initialMass
      ) {
        aboveDouble += 1;
      }

      if (
        Number.isFinite(node.data.steadyMass) &&
        node.data.hasMeaningfulSteadyState &&
        node.data.steadyMass > 0 &&
        node.data.currentMass >= 0.5 * node.data.steadyMass
      ) {
        aboveHalfSteady += 1;
      }
    }

    return { totalMass, steadyStateMass, highestRatio, aboveDouble, aboveHalfSteady };
  }, [nodes]);

  const formattedCurrentTime = useMemo(() => {
    return formatTimeValue(currentTimeValue, timeUnit);
  }, [currentTimeValue, timeUnit]);

  const formattedMaxTime = useMemo(() => {
    return formatTimeValue(maxTimeValue, timeUnit);
  }, [maxTimeValue, timeUnit]);

  const plantDoublingTime = Number(surrogatePredictions.plantDoublingTime?.value);
  const doublingTimeIndex = useMemo(() => {
    return findClosestTimeIndex(timeSeries, plantDoublingTime);
  }, [plantDoublingTime, timeSeries]);
  const hasDoublingTime = Number.isFinite(plantDoublingTime) && doublingTimeIndex >= 0;
  const doublingTimeProgress =
    hasDoublingTime && maxTime > 0 ? Math.max(0, Math.min(100, (doublingTimeIndex / maxTime) * 100)) : 0;
  const doublingTimeValue = hasDoublingTime ? timeSeries[doublingTimeIndex] : Number.NaN;
  const isNearDoublingTime =
    hasDoublingTime && Math.abs(boundedTimeIndex - doublingTimeIndex) <= Math.max(2, Math.ceil(maxTime * 0.008));
  const hasReachedDoublingTime = hasDoublingTime && boundedTimeIndex >= doublingTimeIndex;
  const doublingAlertText = hasDoublingTime
    ? hasReachedDoublingTime
      ? `Plant doubling time reached at ${formatTimeValue(doublingTimeValue, timeUnit)}`
      : `Plant doubling time marker: ${formatTimeValue(doublingTimeValue, timeUnit)}`
    : '';

  const exportableLayout = useMemo(() => {
    const savedPositions = baseNodes
      .filter((node) => node.id !== 'wall_loss')
      .map((node) => ({
        id: node.id,
        position: {
          x: Number(node.position.x),
          y: Number(node.position.y),
        },
      }));

    const wallLossNode = baseNodes.find((node) => node.id === 'wall_loss');

    return {
      savedPositions,
      wallLossPosition: wallLossNode
        ? {
            x: Number(wallLossNode.position.x),
            y: Number(wallLossNode.position.y),
          }
        : null,
    };
  }, [baseNodes]);

  const copyLayoutToClipboard = useCallback(async () => {
    const payload = JSON.stringify(exportableLayout, null, 2);

    try {
      await navigator.clipboard.writeText(payload);
      setLayoutStatus('Layout and wall-loss position copied to clipboard.');
    } catch {
      setLayoutStatus('Clipboard copy failed. Use the text output from the browser console.');
      console.log(payload);
    }

    window.setTimeout(() => {
      setLayoutStatus('');
    }, 2500);
  }, [exportableLayout]);

  const edges = useMemo(() => {
    return allEdges.filter((edge) => {
      if (edge.data?.isLoss && !showLossEdges) {
        return false;
      }

      if (edge.data?.isUnknownFlow && !showUnspecifiedEdges) {
        return false;
      }

      return true;
    });
  }, [allEdges, showLossEdges, showUnspecifiedEdges]);

  return (
    <div
      className="app-shell"
      style={{
        backgroundColor: CANVAS_BASE,
      }}
    >
      <ReactFlow
        nodes={nodes}
        edges={edges}
        nodeTypes={nodeTypes}
        onNodesChange={onNodesChange}
        fitView
        fitViewOptions={{ padding: 0.16 }}
        minZoom={0.2}
        maxZoom={1.6}
        defaultEdgeOptions={{ zIndex: 0 }}
      >
        <Background
          variant={BackgroundVariant.Lines}
          gap={32}
          size={1}
          color={CANVAS_GRID}
        />
        <MiniMap
          pannable
          zoomable
          maskColor="rgba(0, 0, 0, 0.24)"
          nodeColor={(node) => (node.data?.isLossNode ? LOSS_NODE_COLOR : FLOW_COLOR)}
          style={{
            background: 'rgba(239, 242, 238, 0.95)',
            border: `1px solid ${SURFACE_BORDER}`,
            borderRadius: 8,
          }}
        />
        <Controls
          style={{
            background: SURFACE,
            color: PANEL_TEXT,
            border: `1px solid ${SURFACE_BORDER}`,
            borderRadius: 8,
            overflow: 'hidden',
          }}
        />

        <Panel position="top-left">
          <div className="dashboard-panel">
            <div className="panel-kicker">
              <span>Fuel-cycle replay</span>
              <span className={`source-pill is-${dataStatus.state}`}>{dataStatus.source}</span>
            </div>
            <div className="panel-title">
              Tritium Fuel Cycle
            </div>

            <div className="source-strip">
              <div>
                <div className="source-label">Data</div>
                <div className="source-message">{dataStatus.message}</div>
              </div>
              {dataStatus.remoteBaseUrl && (
                <div className="source-url">{dataStatus.remoteBaseUrl}</div>
              )}
              <button
                type="button"
                onClick={() => loadVisualizationData(true)}
                className="fallback-button"
                disabled={dataStatus.state === 'loading' || dataStatus.source === 'Bundled override'}
              >
                {dataStatus.source === 'Bundled override'
                  ? 'Using bundled fallback'
                  : dataStatus.source.startsWith('Remote')
                    ? 'Switch to bundled fallback'
                    : 'Use bundled fallback'}
              </button>
            </div>

            <div className="transport-row">
              <button
                type="button"
                onClick={() => setIsPlaying((playing) => !playing)}
                className="control-button is-primary"
              >
                {isPlaying ? 'Pause' : 'Play'}
              </button>
              <button
                type="button"
                onClick={() => {
                  setIsPlaying(false);
                  setTimeIndex(0);
                }}
                className="control-button"
              >
                Reset
              </button>
              <button
                type="button"
                onClick={copyLayoutToClipboard}
                className="control-button"
              >
                Copy layout
              </button>
            </div>

            {layoutStatus && (
              <div className="layout-status">{layoutStatus}</div>
            )}

            <div className="timeline-block">
              <div className="timeline-labels">
                <span>{formattedCurrentTime}</span>
                <span>{formattedMaxTime} max</span>
              </div>
              <div className="timeline-shell">
                <input
                  className="timeline-range"
                  type="range"
                  min={0}
                  max={maxTime}
                  value={boundedTimeIndex}
                  onChange={(event) => {
                    const nextValue = Number(event.target.value);
                    startTransition(() => setTimeIndex(nextValue));
                  }}
                />
                {hasDoublingTime && (
                  <button
                    type="button"
                    className={`doubling-marker ${hasReachedDoublingTime ? 'is-reached' : ''} ${isNearDoublingTime ? 'is-near' : ''}`}
                    style={{ '--doubling-left': `${doublingTimeProgress}%` }}
                    onClick={() => {
                      setIsPlaying(false);
                      setTimeIndex(doublingTimeIndex);
                    }}
                    aria-label="Jump to plant doubling time"
                    title={formatTimeValue(doublingTimeValue, timeUnit)}
                  >
                    <span />
                  </button>
                )}
              </div>
              <div className="timeline-index">
                sample {boundedTimeIndex.toLocaleString()} / {maxTime.toLocaleString()}
              </div>
              {hasDoublingTime && (
                <div className={`doubling-alert ${hasReachedDoublingTime ? 'is-reached' : ''} ${isNearDoublingTime ? 'is-near' : ''}`}>
                  <span className="doubling-alert-dot" />
                  <span>{doublingAlertText}</span>
                </div>
              )}
            </div>

            <div className="segmented-control">
              {PLAYBACK_OPTIONS.map((option, index) => (
                <button
                  key={option.label}
                  type="button"
                  onClick={() => setPlaybackMode(index)}
                  className={index === playbackMode ? 'is-active' : ''}
                >
                  {option.label}
                </button>
              ))}
            </div>

            <div className="metric-grid">
              <MetricCard label="Total current mass" value={formatMassWithUnit(massStats.totalMass)} accent={FLOW_COLOR} />
              <MetricCard label="Total steady-state mass" value={formatMassWithUnit(massStats.steadyStateMass)} accent="#d9e0e5" />
              <MetricCard label="Above 2x init" value={String(massStats.aboveDouble)} accent={WARNING_RED} />
              <MetricCard label="Above 0.5 ss" value={String(massStats.aboveHalfSteady)} accent={WARNING_AMBER} />
            </div>

            <div className="panel-section surrogate-section">
              <div className="panel-section-title">
                Surrogate predictions
              </div>
              <div className="surrogate-grid">
                {Object.entries(surrogatePredictions).map(([key, prediction]) => (
                  <PredictionCard
                    key={key}
                    label={prediction.label}
                    value={formatPredictionValue(prediction)}
                  />
                ))}
              </div>
            </div>

            <div className="panel-section">
              <div className="panel-section-title">
                Edge layers
              </div>
              <label className="toggle-row">
                <input
                  type="checkbox"
                  checked={showLossEdges}
                  onChange={(event) => setShowLossEdges(event.target.checked)}
                />
                <span>Show loss edges</span>
              </label>
              <label className="toggle-row">
                <input
                  type="checkbox"
                  checked={showUnspecifiedEdges}
                  onChange={(event) => setShowUnspecifiedEdges(event.target.checked)}
                />
                <span>Show unweighted injector links</span>
              </label>
            </div>

            <div className="panel-section">
              <div className="panel-section-title">
                Edge legend
              </div>
              <LegendRow color={FLOW_COLOR} label="Process line. Thicker stroke indicates larger fractional inflow." />
              <LegendRow color={UNSPECIFIED_COLOR} label="Injector link exists, but no fractional split was exported." />
              <LegendRow color={LOSS_EDGE_COLOR} label="Wall-loss branch routed to the sink node." />
            </div>
            <div className="threshold-key">
              <ThresholdKeyItem color={WARNING_AMBER} label="0.5 ss" />
              <ThresholdKeyItem color={WARNING_RED} label="2x init" />
              <ThresholdKeyItem color={WARNING_PURPLE} label="both" />
            </div>
          </div>
        </Panel>
      </ReactFlow>
    </div>
  );
}

function MetricCard({ accent, label, value }) {
  return (
    <div
      className="metric-card"
      style={{
        '--metric-accent': accent,
      }}
    >
      <div className="metric-label">{label}</div>
      <div className="metric-value">{value}</div>
    </div>
  );
}

function PredictionCard({ label, value }) {
  return (
    <div className="prediction-card">
      <div className="prediction-label">{label}</div>
      <div className="prediction-value">{value}</div>
    </div>
  );
}

function ThresholdKeyItem({ color, label }) {
  return (
    <span className="threshold-key-item" style={{ '--threshold-color': color }}>
      <span className="threshold-key-dot" />
      {label}
    </span>
  );
}

function LegendRow({ color, dashed = false, label }) {
  return (
    <div className="legend-row">
      <span
        className="legend-line"
        style={{
          '--legend-color': color,
          '--legend-style': dashed ? 'dashed' : 'solid',
        }}
      />
      <span>{label}</span>
    </div>
  );
}
