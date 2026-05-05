from __future__ import annotations

import argparse
import asyncio
import os
import sys

from fastmcp import Client


DEFAULT_MCP_URL = os.environ.get("RHINO_MCP_URL", "http://127.0.0.1:8000/mcp")


async def run_demo(args: argparse.Namespace) -> None:
    client = Client(args.mcp_url)

    async with client:
        result = await client.call_tool(
            "build_graph_for_nearest_simulation",
            {
                "Ndotminus": args.Ndotminus,
                "beta": args.beta,
                "write_visualization_files": args.write_json_files,
                "start_visualization": not args.no_start_visualization,
                "visualization_port": args.visualization_port,
            },
        )

    data = result.data
    nearest = data["nearest_simulation"]
    visualization = data["visualization"]
    graph = data["graph"]

    print("RHINO surrogate demo")
    print("====================")
    print(f"Inputs:")
    print(f"  Ndotminus: {args.Ndotminus} g/d")
    print(f"  beta:      {args.beta}")
    print()
    print("Nearest simulation:")
    print(f"  id:       {nearest['simulation_id']}")
    print(f"  path:     {nearest['path']}")
    print(f"  exists:   {nearest['path_exists']}")
    print(f"  distance: {nearest['distance_normalized']:.6g}")
    print()
    print("Graph:")
    print(f"  nodes: {len(graph['nodes'])}")
    print(f"  edges: {len(graph['edges'])}")
    print()
    print("Visualization:")
    print(f"  status: {visualization['status']}")
    if visualization.get("url"):
        print(f"  url:    {visualization['url']}")
    if visualization.get("reason"):
        print(f"  reason: {visualization['reason']}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the RHINO demo from two physical inputs to visualization."
    )
    parser.add_argument(
        "--Ndotminus",
        type=float,
        required=True,
        help="Tritium burning rate in grams per day.",
    )
    parser.add_argument(
        "--beta",
        type=float,
        required=True,
        help="Burn fraction as a unitless fraction.",
    )
    parser.add_argument(
        "--mcp-url",
        default=DEFAULT_MCP_URL,
        help=f"MCP server URL. Defaults to {DEFAULT_MCP_URL!r}.",
    )
    parser.add_argument(
        "--visualization-port",
        type=int,
        default=5173,
        help="Local Vite visualization port.",
    )
    parser.add_argument(
        "--no-start-visualization",
        action="store_true",
        help="Build live graph data but do not start the visualization server.",
    )
    parser.add_argument(
        "--write-json-files",
        action="store_true",
        help="Also write graph JSON files into react_flow_viz/public.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    try:
        asyncio.run(run_demo(parse_args()))
    except RuntimeError as exc:
        if "Client failed to connect" not in str(exc):
            raise

        print(
            "Could not connect to the RHINO MCP server.\n\n"
            f"The client tried: {DEFAULT_MCP_URL}\n\n"
            "Start the server in another terminal first:\n"
            "  cd RHINO/surrogate\n"
            "  conda activate scspdemo\n"
            "  python mcp_server_rhino.py\n\n"
            "If the server is using a different port, pass it explicitly, for example:\n"
            "  python demo_client.py --Ndotminus 500 --beta 0.1 "
            "--mcp-url http://127.0.0.1:8001/mcp",
            file=sys.stderr,
        )
        sys.exit(1)
