import asyncio
from fastmcp import Client

client = Client("mcp_server_rhino.py")

async def main():
    async with client:
        tools = await client.list_tools()
        print("Available tools:")
        for tool in tools:
            print("-", tool.name)

        result = await client.call_tool(
            "predict_rhino_surrogate",
            {
                "I0_SD": 100.0,
                "Ndotminus": 500.0,
                "beta": 0.1,
            },
        )
        print("\nTool result:")
        print(result)

        nearest = await client.call_tool(
            "find_nearest_simulation",
            {
                "I0_SD": 100.0,
                "Ndotminus": 500.0,
                "beta": 0.1,
            },
        )
        print("\nNearest simulation:")
        print(nearest.data)


        result_and_nearest = await client.call_tool(
            "predict_and_compare_to_nearest_simulation",
            {
                "I0_SD": 100.0,
                "Ndotminus": 500.0,
                "beta": 0.1,
            },
        )
        
        print("\nStructured Python object:")
        print(result_and_nearest.data)
        
        print("\nRaw structured JSON:")
        print(result_and_nearest.structured_content)
        
        if result_and_nearest.data is None:
            print("\nFallback content blocks:")
            for block in result_and_nearest.content:
                if hasattr(block, "text"):
                    print(block.text)

        result_graph = await client.call_tool(
            "build_graph_for_nearest_simulation",
            {
                "I0_SD": 100.0,
                "Ndotminus": 500.0,
                "beta": 0.1,
            },
        )

        data = result_graph.data
        print(data.keys())

        graph = data["graph"]
        nodes = graph["nodes"]
        edges = graph["edges"]
        time = graph["time"]

        print("n nodes:", len(nodes))
        print("n edges:", len(edges))
        print("time info:", time)        

asyncio.run(main())